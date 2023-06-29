#include "eulerFunc.H"
#include "diffusionFunc.H"
#include "thermoTransport.H"
#include "particleFunc.H"
#include "AmrLevelAdv.H"
// #include "constants.H"
// #include "constants.H"
#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

#include <cmath>
#include <iostream>
#include <algorithm> 

extern int enIC;
extern int NUM_STATE;
extern int particle;
extern int n_cell;
extern const int spacedim;
extern double R, pi, one_atm_Pa;
extern double dp0;               // mu m                 Initial particle diameter
extern double delta0;            // [-]                  Initial oxide layer thickness to particle size ratio
extern double rhoFe;             // kg/m^3               Solid-phase iron(Fe) density
extern double rhoFeO;            // kg/m^3               Solid-phase FeO density
extern double rhoFe3O4;          // kg/m^3               Solid-phase Fe3O4 density
extern double M_Fe;              // kg/mol               Molar mass of iron (Fe)
extern double M_FeO;             // kg/mol               Molar mass of wustite (FeO)
extern double M_Fe3O4;           // kg/mol               Molar mass of wustite (FeO)
extern double M_Fe2O3;           // kg/mol               Molar mass of wustite (FeO)
extern double HfuFe;             // J/kg                 Heat of fusion to melt Fe
extern double HfuFeO;            // J/kg                 Heat of fusion to melt FeO
extern double HfuFe3O4;          // J/kg                 Heat of fusion to melt Fe3O4
extern double k0FeOs;            // m^2/s                Pre-exponential factor for kinetic rate of FeO oxide layer growth
extern double k0Fe3O4s;          // m^2/s                Pre-exponential factor for kinetic rate of Fe3O4 oxide layer growth
extern double TaFeOs;            // K                    Activation temperature for kinetic rate of FeO oxide layer growth
extern double TaFe3O4s;          // K                    Activation temperature for kinetic rate of Fe3O4 oxide layer growth
extern double nFeFeO;
extern double nFeFe3O4;
extern double nO2FeO;
extern double nO2Fe3O4;
extern double TpInitial;
extern double M_O2, M_N2;
extern double qFeOs,qFeOl,qFe3O4s,qFe2O3s,qFeOg;
extern double meltFe,meltFeO,meltFe3O4;
extern double mFe0,mFeO0,mFe3O40,interDist;


using namespace amrex;

// Use this file to write functions required for particle calculations

void
AmrLevelAdv::initParticles (const MultiFab& S_new)
{
    const int lev = 0;
    Real patch = 0;
    int totalParIter = 0;

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    //    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        Real Np = 0;
        const Box& bx = mfi.tilebox();
        const Dim3 lo = lbound(bx);
        const Dim3 hi = ubound(bx);
        const Real* dx = geom.CellSize();
        double dX = dx[0];

        std::ifstream locations;

        if ((enIC==12)||(enIC==13)){
            locations.open("setupScripts/locationSingle.txt", std::ios_base::in);
        }
        else if (enIC==14){
            locations.open("setupScripts/locations.txt", std::ios_base::in);
        }

        auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];

        Real loxBound, hixBound, loyBound, hiyBound;

        loxBound = lo.x;
        hixBound = hi.x;
        
        Vector<double> x_coord;
        Vector<double> y_coord;

        float x_txt,y_txt;

        if (spacedim == 2){
            Real loyBound = lo.y;
            Real hiyBound = hi.y;
            while (locations >> x_txt >> y_txt)
            {
                if ((x_txt < hi.x)&&(x_txt >= lo.x)) {
                    if ((y_txt < hi.y)&&(y_txt >= lo.y)){
                        x_coord.push_back (x_txt*dX);
                        y_coord.push_back (y_txt*dX);
                        Np += 1;
                        // std::cout << "x, y from locations: " << x_txt << " " << y_txt << std::endl;
                    }
                }
            }
        }    
        else {
            while (locations >> x_txt)
            {
                if ((x_txt < hi.x)&&(x_txt >= lo.x)) {
                    x_coord.push_back (x_txt*dX);
                    Np += 1;
                    std::cout << "x from location: " << x_txt << std::endl;
                }
            }
        }

        // std::cout << "total number of particles in this patch: "  << Np << std::endl;
        
        for (int i = 0; i < Np; i++){

            ParticleType p;
            p.id()   = ParticleType::NextID();
            p.cpu()  = ParallelDescriptor::MyProc();
            
            p.pos(0) = x_coord[i];
            // std::cout << "particle position: " << p.pos(0) << std::endl;
            if (spacedim == 2){
                p.pos(1) = y_coord[i];
            }

            double energy0;
            particleInit(energy0);
            
            if (enIC==14){
                if ((p.pos(0)/dX) < n_cell*0.1){
                    energy0 = Hparticle(mFe0,mFeO0,mFe3O40,1270,0,0,0);
                }
            }    

            p.rdata(RealData::u)     = 0.0;         // Store the x-velocity here (unused if simulation is 1-D)
            p.rdata(RealData::v)     = 0.0;       // Store the y-velocity here (unused if simulation is 1-D)
            p.rdata(RealData::w)     = 0.0;       // Store the z-velocity here (unused if simulation is 1-D)
            p.rdata(RealData::mFe)    = mFe0;      // Store the initial Fe mass here
            p.rdata(RealData::mFeO)   = mFeO0;     // Store the initial FeO mass here
            p.rdata(RealData::mFe3O4) = mFe3O40;   // Store the initial Fe3O4 mass here
            p.rdata(RealData::Hp)     = energy0;   // Store the initial particle enthalpy here
            p.rdata(RealData::LFe)    = 0.0;       // Store the Fe melt progress variable here
            p.rdata(RealData::LFeO)   = 0.0;       // Store the FeO melt progress variable here
            p.rdata(RealData::LFe3O4) = 0.0;       // Store the Fe3O4 melt progress variable here

            p.idata(IntData::Fe)     = 0;         // Store the Fe melt flag variable here
            p.idata(IntData::FeO)    = 0;         // Store the FeO melt flag variable here
            p.idata(IntData::Fe3O4)  = 0;         // Store the Fe3O4 melt flag variable here
            p.idata(IntData::pIter)  = totalParIter; // Store the particle counter here
            p.idata(IntData::regime) = 0;            // Store the particle combustion regime

            totalParIter += 1;

            particles.push_back(p);

            // std::cout << "x, y: " << p.pos(0) << " " << p.pos(1) << std::endl;

        }
        // std::cout << "patch number: " << patch << std::endl;
    }
    Redistribute();
}

Vector<double> 
AmrLevelAdv::getParticleInfo(Vector<double>& pReal, Vector<int>& pInt)
{
  const int lev = 0;
  Real x, y=0, up, vp, wp, mFe, mFeO, mFe3O4, Hp, Tp, LFe, LFeO, LFe3O4;
  int  phaseFe, phaseFeO, phaseFe3O4;
  Vector<double> data(3);

//   std::cout << "in getParticleInfo" << std::endl;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){
    // for (MFIter mfi(S_new); mfi.isValid(); ++mfi)    
    // {
  
    const int grid_id = mfi.index();
    const int tile_id = mfi.LocalTileIndex();
    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    auto& particles = particle_tile.GetArrayOfStructs();
    
    const int np = particles.numParticles();
    // std::cout << "number of particles in grid in getParticleInfo, within advance: " << np << std::endl;

    for(int pindex = 0; pindex < np; ++pindex) {
        ParticleType& p = particles[pindex];
        const IntVect& iv = this->Index(p, lev);
        x = p.pos(0);
        if (spacedim == 2){
            y = p.pos(1);
        }
        for (int h = 0; h < RealData::ncomps; h++){
            pReal[h] = p.rdata(h);
        }
        for (int h = 0; h < IntData::ncomps; h++){
            pInt[h] = p.idata(h);
        }
        up     = pReal[RealData::u];
        vp     = pReal[RealData::v];
        wp     = pReal[RealData::w];
        mFe    = pReal[RealData::mFe];
        mFeO   = pReal[RealData::mFeO];
        mFe3O4 = pReal[RealData::mFe3O4];
        Hp     = pReal[RealData::Hp];
        LFe    = pReal[RealData::LFe];
        LFeO   = pReal[RealData::LFeO];
        LFe3O4 = pReal[RealData::LFe3O4];
        phaseFe    = pInt[IntData::Fe];
        phaseFeO   = pInt[IntData::FeO];
        phaseFe3O4 = pInt[IntData::Fe3O4];
        Tp = Tparticle(mFe,mFeO,mFe3O4,Hp,phaseFe,phaseFeO,phaseFe3O4,LFe,LFeO,LFe3O4);
    }
  }
  data[0] = x;
  data[1] = y;
  data[2] = Tp;
  return data;
}

void 
AmrLevelAdv::updateParticleInfo(MultiFab& Sborder, const double& dt, const double& dx, const double& dy)
{
  const int lev = 0;
  Real x, y, lx, ly, vCell;
  int  i=0, j=0, k=0;
  Vector<double> q(NUM_STATE,0);
  Vector<double> qSource(NUM_STATE,0);
  Vector<double> pSource(RealData::ncomps,0);
  Vector<double> pReal(RealData::ncomps,0);
  Vector<int>    pInt(IntData::ncomps,0);

  
//   std::cout << "in updateParticleInfo" << std::endl;

  for (MFIter mfi(Sborder); mfi.isValid(); ++mfi){

    const int grid_id = mfi.index();
    const int tile_id = mfi.LocalTileIndex();
    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    auto& particles = particle_tile.GetArrayOfStructs();
    const int np = particles.numParticles();
    // std::cout << "number of particles in grid in updateParticleInfo, within advance: " << np << std::endl;

    // Access S_border multifab information here
    const Box& bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);
    const auto& arr = Sborder.array(mfi);

    for(int pindex = 0; pindex < np; ++pindex) {
        ParticleType& p = particles[pindex];
        const IntVect& iv = this->Index(p, lev);

        x  = p.pos(0);   //given in meters, lo.x is in cell count, so scale by dx to get meters
        i  = static_cast<int>(Math::floor(x/dx));
        // std::cout << "x-position: " << x << std::endl;
        // std::cout << "cell number is: " << i << std::endl;
        // std::cout << "lo.x and hi.x are: " << lo.x << ", " << hi.x << std::endl; 

        if (spacedim == 2){
            y  = p.pos(1);
            j  = static_cast<int>(Math::floor(y/dx));
            // std::cout << "y-position: " << y << std::endl;
            // std::cout << "position in y cell number is: " << y/dy << ", cell number is: " << j << std::endl;
        }
        
        // std::cout << "rho: " << arr(i,j,k,0) << std::endl;
        // std::cout << "x- and y-momentum: " << arr(i,j,k,1) << " " << arr(i,j,k,2) << std::endl;
        // std::cout << "energy: " << arr(i,j,k,3) << std::endl;
        // std::cout << "O2 and N2 concentration: " << arr(i,j,k,4) << " " << arr(i,j,k,5) << std::endl;

        for (int h = 0; h < NUM_STATE; h++){
            q[h] = arr(i,j,k,h);
            if (q[h] != q[h]){
                std::cout << "Nan found in updateParticle for gas-phase, variable h: " << h << std::endl;
                std::cout << "i,j,k: " << i << " " << j << " " << k << std::endl;
                std::cout << "grid_id: " << grid_id << ", tile_id: " << tile_id << std::endl;
                Abort("nan found");
            }
        }
        for (int h = 0; h < RealData::ncomps; h++){
            pReal[h] = p.rdata(h);
        }
        for (int h = 0; h < IntData::ncomps; h++){
            pInt[h] = p.idata(h);
        }

        getSource(qSource,pSource,q,pReal,pInt,dt,dx,dy);

        p.rdata(RealData::Hp)     = pReal[RealData::Hp];
        p.idata(IntData::Fe)      = pInt[IntData::Fe];
        p.idata(IntData::FeO)     = pInt[IntData::FeO];
        p.idata(IntData::Fe3O4)   = pInt[IntData::Fe3O4];
        p.idata(IntData::regime)   = pInt[IntData::regime];
        
        p.pos(0) = p.pos(0) + dt*p.rdata(RealData::u);
        if (spacedim == 2){
            p.pos(1) = p.pos(1) + dt*p.rdata(RealData::v);
        }
        for (int h = 0; h < RealData::ncomps; h++){
            p.rdata(h) = p.rdata(h) + dt*pSource[h];
        }

        if (particle == 2){ // if we enable two-way coupling
            vCell = dx*interDist*interDist;
            for (int h = 0; h < NUM_STATE; h++){
                    
                if (qSource[h] != qSource[h]){
                    std::cout << "Nan BEFORE applying LagSource in qSource, \nrho rhou rhov e o2 n2" << arr(i,j,k,0) << " " << arr(i,j,k,1) << \
                    " " << arr(i,j,k,2) << " " << arr(i,j,k,3) << " " << arr(i,j,k,4) << " " << arr(i,j,k,5) << std::endl;
                    std::cout << "qSrc: rho rhou rhov e o2 n2" << qSource[0] << " " << qSource[1] << \
                    " " << qSource[2] << " " << qSource[3] << " " << qSource[4] << " " << qSource[5] << std::endl;
                    Abort("nan found before applying lagrangian source");
                }
                
                q[h] = q[h] + dt*(1.0/vCell)*qSource[h];
                arr(i,j,k,h) = q[h];

                if (arr(i,j,k,h) != arr(i,j,k,h)){
                    std::cout << "cell number: " << i << std::endl;
                    std::cout << "dx: " << dx << ", interDist: " << interDist << std::endl;
                    std::cout << "dt: " << dt << ", vCell: " << vCell << std::endl;
                    std::cout << "Gas variables before update, \nrho rhou rhov e o2 n2" << q[0] << " " << q[1] << \
                    " " << q[2] << " " << q[3] << " " << q[4] << " " << q[5] << std::endl;
                    std::cout << "Nan AFTER applying LagSource, \nrho rhou rhov e o2 n2" << arr(i,j,k,0) << " " << arr(i,j,k,1) << \
                    " " << arr(i,j,k,2) << " " << arr(i,j,k,3) << " " << arr(i,j,k,4) << " " << arr(i,j,k,5) << std::endl;
                    std::cout << "qSrc: rho rhou rhov e o2 n2\n" << qSource[0] << " " << qSource[1] << \
                    " " << qSource[2] << " " << qSource[3] << " " << qSource[4] << " " << qSource[5] << std::endl;
                    std::cout << "Nan found in particle to gas update, variable h: " << h << std::endl;
                    Abort("nan found after particle to gas update");
                }
            }
        }
        
        
    }
    // Redistribute();
  }
  Redistribute();
}

void getSource(Vector<double>& qSource, Vector<double>& pSource, const Vector<double>& q, Vector<double>& pReal, \
               Vector<int>& pInt, const double& dt, const double& dx, const double& dy){
    
    // Declare particle and gas variables
    Real up, vp, wp, mFe, mFeO, mFe3O4, Hp, LFe, LFeO, LFe3O4;
    int  phaseFe, phaseFeO, phaseFe3O4;
    Real rho, rhou, rhov, energy, rhoYO2, rhoYN2, u, v, p, Tgas, YO2, YN2;

    // Declare variables used for vapor pressure calculation
    Real pFe=0, pFeO=0, pO2=0, pN2=0, XO2=0, XN2=0, XFep=0, XFeOp=0, XO2p=0, XN2p=0, YFep=0, YFeOp=0, YO2p=0, YN2p=0;
    Vector<double> yGas(gases::ncomps);

    // Declare variables used for drag calculation
    Real Tp=0, vTot=0, rp=0, dp=0, rhop=0, Re=0, mu=0, uRel=0, vRel=0, CD=0;
    Real rhoFel=rhoFe, rhoFeOl=rhoFeO;

    // Declare variables used for evaporation calculations
    Real mdotFe=0, mdotFeO=0, mdotFeOevap=0, qFeO;

    // Declare variables used for heat transfer calculations
    Real Nu=0, Ap=0, epsilon=0, SB=0, kgas=0, conv=0, rad=0, Pr=0, NuSt=0, evap=0;

    // Declare variables for external transport rates
    Real DO2=0, DFe=0, DFeO=0, rhoO2g=0, rhoFeg=0, rhoFeOg=0, ScO2=0, ScFe=0, ScFeO=0, ShO2=0, ShFe=0, ShFeO=0;
    Real mdotO2=0, mdotO2k=0, mdotO2d=0, ShSt=0, cpgas=0, phi=0, Bm=0, Bt=0, Le=0, psi=0;
    Vector<double> mixDiffCoeffs(gases::ncomps);

    // Declare boundary-layer averaged quantities
    Real pfilm=0, Tfilm=0, YO2film=0, YN2film=0, YFefilm=0, YFeOfilm=0, rhofilm=0, Mfilm=0;

    // Declare source terms
    Real dupdt=0, dvpdt=0, dHpdt=0, dmO2dt=0, dmFeOformdt=0, dmFe3O4formdt=0, dHgOxidt=0; //, dHpdt, dmFedt, dmFeOdt, dmFe3O4dt;
    Vector<double> dmdt(3);

    rho    = q[gasVar::rho];
    rhou   = q[gasVar::rhou];
    rhov   = q[gasVar::rhov];
    energy = q[gasVar::E];
    rhoYO2 = q[gasVar::rhoYO2];
    rhoYN2 = q[gasVar::rhoYN2];

    // std::cout << "rho rhou ener O2 N2: " << rho << " " << rhou << " " << energy << " " << rhoYO2/rho << " " << rhoYN2/rho << std::endl;

    u    = rhou/rho;
    v    = rhov/rho;
    YO2  = rhoYO2/rho;
    YN2  = rhoYN2/rho;
    Tgas = Tg(rho,u,v,YO2,YN2,energy);
    p    = pressure(rho,YO2,YN2,Tgas);

    // std::cout << "Tgas, p: " << Tgas << " " << p << std::endl;

    up     = pReal[RealData::u];
    vp     = pReal[RealData::v];
    wp     = pReal[RealData::w];
    mFe    = pReal[RealData::mFe];
    mFeO   = pReal[RealData::mFeO];
    mFe3O4 = pReal[RealData::mFe3O4];
    Hp     = pReal[RealData::Hp];
    LFe    = pReal[RealData::LFe];
    LFeO   = pReal[RealData::LFeO];
    LFe3O4 = pReal[RealData::LFe3O4];
    phaseFe    = pInt[IntData::Fe];
    phaseFeO   = pInt[IntData::FeO];
    phaseFe3O4 = pInt[IntData::Fe3O4];
    Tp         = Tparticle(mFe,mFeO,mFe3O4,Hp,phaseFe,phaseFeO,phaseFe3O4,LFe,LFeO,LFe3O4);
    Tfilm      = filmAverage(Tp,Tgas);

    // std::cout << "mFe, mFeO: " << mFe << " " << mFeO << "\n" << std::endl;
    // std::cout << "mFe3O4, Hp: " << mFe3O4 << " " << Hp << "\n" << std::endl;
    // std::cout << "Tp, Tg: " << Tp << " " << Tgas << "\n" << std::endl;

    // With the particle temperature known, we first compute the vapor pressures of gas-phase Fe and FeO
    // resulting from the gas-liquid equilibrium at the particle surface. If the particle temperature is 
    // below the melting point of Fe3O4, the vapor pressure is zero (use initialized value).
    if (Tp > 1870){
        pFe  = one_atm_Pa*pow(10.0,6.041-2.095*1e4/Tp);  // vapor pressure of Fe evaporating from liquid FeO
        pFeO = one_atm_Pa*pow(10.0,5.962-2.175*1e4/Tp);  // vapor pressure of FeO evaporating from liquid FeO
        pO2  = one_atm_Pa*pow(10.0,5.361-2.016*1e4/Tp);  // vapor pressure of O2 released from liquid FeO
        rhoFeOl = 4.35*1.0e3;
        rhoFel  = 1e3*(8.523-(8.358e-4)*Tp);
    }
    XO2 = (YO2/M_O2)/(YO2/M_O2+YN2/M_N2);  // mole fraction of O2 in gas-phase 
    XN2 = (YN2/M_N2)/(YO2/M_O2+YN2/M_N2);  // mole fraction of N2 in gas-phase 
    // pN2 = one_atm_Pa-pFe-pFeO-pO2;
    pN2 = XN2*p;                           // partial pressure of N2 at the particle surface

    XO2p = pO2/(pO2+pN2+pFe+pFeO);
    XN2p = pN2/(pO2+pN2+pFe+pFeO);
    XFep = pFe/(pO2+pN2+pFe+pFeO);
    XFeOp= pFeO/(pO2+pN2+pFe+pFeO);

    // std::cout << "o2 n2 Xp: " << XO2p << " " << XN2p << std::endl;

    pfilm = filmAverage(pO2+pN2+pFe+pFeO,p);

    // std::cout << "film pressure: " << pfilm << std::endl;

    yGas = getMassFractions(XO2p, XN2p, XFep, XFeOp);
    YO2p = yGas[gases::O2];
    YN2p = yGas[gases::N2];
    YFep = yGas[gases::Fe];
    YFeOp= yGas[gases::FeO];

    // std::cout << "o2 n2 Yp: " << YO2p << " " << YN2p << std::endl;

    // Calculate boundary-layer averaged mass fractions and density:
    YO2film = filmAverage(YO2p,YO2);
    YN2film = filmAverage(YN2p,YN2);
    YFefilm = filmAverage(YFep,0);
    YFeOfilm = filmAverage(YFeOp,0);
    Mfilm    = 1.0/(YO2film/M_O2+YN2film/M_N2+YFefilm/M_Fe+YFeOfilm/M_FeO);
    rhofilm  = pfilm*Mfilm/(R*Tfilm);

    // std::cout << "rhofilm " << rhofilm << std::endl;

    // Calculate necessary particle parameters here:
    vTot = mFe/rhoFel + mFeO/rhoFeOl + mFe3O4/rhoFe3O4;    // total particle volume, m^3
    rp   = pow(3.0*vTot*0.25/pi,1.0/3.0);                // particle outer radius, m
    dp   = 2*rp;                                     
    rhop = (mFe+mFeO+mFe3O4)/vTot;                       // average particle density, kg/m^3

    // std::cout << "particle radius: " << rp << std::endl;
    
    // Calculate Reynolds number for particle flow here:
    uRel = fabs(u-up);                               // relative velocity (slip velocity), m/s
    mu = muMix(muO2(Tfilm),muN2(Tfilm),muFe(Tfilm),muFeO(Tfilm),YO2film,YN2film,YFefilm,YFeOfilm);  // mixture viscosity of the gas
    Re = 2.0*rp*rhofilm*uRel/mu;                         // Reynolds number of the particle flow

    // std::cout << "mugas: " << mu << std::endl;

    // Calculate drag coefficient here:
    CD = Cdrag(Re);

    // Calculate the rate of change of particle velocity here:
    dupdt = (3*CD*rhofilm/(8*rp*rhop))*(u-up)*uRel; // + 9.81;
    pSource[RealData::u] = dupdt; 

    if (spacedim == 2){ // process repeated in y-direction for a 2-D simulation
        vRel  = fabs(v-vp);
        Re    = 2.0*rp*rhofilm*vRel/mu;
        CD    = Cdrag(Re);
        dvpdt = (3*CD*rhofilm/(8*rp*rhop))*(v-vp)*vRel;
        pSource[RealData::v] = dvpdt;
    }

    // Calculate mass diffusivity and oxygen density here:
    mixDiffCoeffs = getMixDiffCoeffs(Tfilm,pfilm,YO2film,YN2film,YFefilm,YFeOfilm); // use film temperature and pressure
    DO2   = mixDiffCoeffs[gases::O2];
    DFe   = mixDiffCoeffs[gases::Fe];
    DFeO  = mixDiffCoeffs[gases::FeO];
    rhoO2g = (pfilm*M_O2)/(R*Tfilm);                    // partial density of O2 in the boundary-layer film, kg/m^3
    rhoFeg = (pfilm*M_Fe)/(R*Tfilm);
    rhoFeOg = (pfilm*M_FeO)/(R*Tfilm);

    // std::cout << "DO2 Fe FeO: " << DO2 << " " << DFe << " " << DFeO << std::endl;
    // std::cout << "pO2 Fe FeO: " << pO2 << " " << pFe << " " << pFeO << std::endl;

    // Calculate Schmidt number and Sherwood number from Froessling equation:
    ScO2  = mu/(rhofilm*DO2);
    ScFe  = mu/(rhofilm*DFe);
    ScFeO = mu/(rhofilm*DFeO);
    phi   = 1.0-1.1*(1.0-pow((Tp/Tgas),0.2));
    ShO2  = 2.0+(0.02*pow(ScO2,1.0/3.0)*sqrt(Re)+0.33*pow(ScO2,1.0/3.0)*pow(Re,2.0/3.0))*phi;
    ShFe  = 2.0+(0.02*pow(ScFe,1.0/3.0)*sqrt(Re)+0.33*pow(ScFe,1.0/3.0)*pow(Re,2.0/3.0))*phi;
    ShFeO = 2.0+(0.02*pow(ScFeO,1.0/3.0)*sqrt(Re)+0.33*pow(ScFeO,1.0/3.0)*pow(Re,2.0/3.0))*phi;

    Bm = (YO2-YO2p)/(YO2p-1.0);
    cpgas = cpMixFe(cpO2(Tfilm),cpN2(Tfilm),cpFeG(Tfilm),cpFeOG(Tfilm),YO2film,YN2film,YFefilm,YFeOfilm);
    kgas  = kMix(kO2(Tfilm),kN2(Tfilm),kFe(Tfilm),kFeO(Tfilm),\
                 YO2film,YN2film,YFefilm,YFeOfilm);
    Pr = mu*cpgas/kgas; 
    Nu = 2.0+(0.02*pow(Pr,1.0/3.0)*sqrt(Re)+0.33*pow(Pr,1.0/3.0)*pow(Re,2.0/3.0))*phi;
    Le = ScO2/Pr;

    psi = (ShO2/Nu)*(1.0/Le)*((cpO2(Tfilm)/M_O2)/cpgas);
    Bt  = pow((1.0+Bm),psi) - 1.0;
    ShSt = ShO2*log(1.0+Bm)/Bm;

    // Calculate external transport rates
    mdotO2  = std::max(2*pi*rp*ShSt*rhoO2g*DO2*(XO2-XO2p),0.0);               // mass transport rate of O2 from the bulk gas to the particle, kg/s
    mdotFe  = 2*pi*rp*ShFe*rhoFeg*DFe*XFep;
    mdotFeO = 2*pi*rp*ShFeO*rhoFeOg*DFeO*XFeOp;
    
    mdotFeOevap = M_FeO*(mdotFeO/M_FeO + mdotFe/M_Fe);
    mdotO2d = mdotO2 - (mdotFe*(3.0/4.0)*M_O2/M_Fe) - (mdotFeO*(1.0/4.0)*M_O2/M_FeO);

    // std::cout << "kmix: " << kgas << std::endl;
    // std::cout << "Re ScO2 rp ShO2 ShSt XO2 XO2p: " << Re << " " << ScO2 << " "  << rp << " " << ShO2 << " " << ShSt << " " << XO2 << std::endl;
    // std::cout << "mdot O2 Fe FeO: " << mdotO2 << " " << mdotFe << " " << mdotFeO << std::endl;
    // std::cout << "mdotO2d: " << mdotO2d << std::endl;

    // Calculate change of Fe, FeO, and Fe3O4 mass based on temperature and oxide layer.
    // Compare the rate at which the kinetics predict a total consumption of oxygen, compare to the molecular diffusion
    // rate. Slower rate is limiting for the overall oxidation reaction.
    dmdt = getOxidationRates(mFe,mFeO,mFe3O4,Tp,rp);
    mdotO2k = dmdt[1]*nO2FeO + dmdt[2]*nO2Fe3O4;     // total kinetic rate of O2 consumption, kg/s

    // std::cout << "kinetic rate: " << mdotO2k << ", diffusion rate: " << mdotO2d << std::endl;
    // std::cout << "fraction of Fe mass remaining: " << mFe/mFe0 << std::endl;
    // std::cout << "fraction of o2 mass fraction remaining: " << YO2/(0.2329) << std::endl;

    // To prevent unboundedness, we only allow the particles to consume oxygen and release heat if
    // there is more than 1% of the initial Fe mass, and more than 1% of the ambient oxygen mole fraction
    // in the gas cell.

    pInt[IntData::regime] = 0;
    if ((mFe/mFe0 > 0.01)&&(YO2/(0.2329) > 0.01)){ 
        if (mdotO2d > mdotO2k){ // if the molecular diffusion rate is FASTER than the kinetic rate of O2 consumption
            // reaction is kinetically-controlled
            pSource[RealData::mFe] = dmdt[0];
            pSource[RealData::mFeO] = dmdt[1];
            pSource[RealData::mFe3O4] = dmdt[2];
            // std::cout << "kinetics-controlled" << std::endl;
            epsilon = 0.88;
            dmO2dt = -nO2FeO*pSource[RealData::mFeO]-nO2Fe3O4*pSource[RealData::mFe3O4];
            dmFeOformdt   = pSource[RealData::mFeO];
            dmFe3O4formdt = pSource[RealData::mFe3O4];
        }
        else{ // if the kinetic rate of O2 consumption is FASTER than the molecular diffusion rate
            // the particle can only consume the O2 supplied from the gas, hence, diffusion-controlled
            if (Tp < 1870){ // if particle is still at least partially in solid-phase, partition the O2 as in Mi et al. (CnF, 2022)
                pSource[RealData::mFeO]   = (mdotO2*dmdt[1]*nO2FeO/(dmdt[1]*nO2FeO+dmdt[2]*nO2Fe3O4))/nO2FeO;
                pSource[RealData::mFe3O4] = (mdotO2*dmdt[2]*nO2Fe3O4/(dmdt[1]*nO2FeO+dmdt[2]*nO2Fe3O4))/nO2Fe3O4;
                epsilon = 0.88;
                dmO2dt = -nO2FeO*pSource[RealData::mFeO]-nO2Fe3O4*pSource[RealData::mFe3O4];
                Nu = Nu*log(1.0+Bt)/Bt;
                dmFeOformdt   = pSource[RealData::mFeO];
                dmFe3O4formdt = pSource[RealData::mFe3O4];
            }
            else{ // if particle is completely in liquid-phase, use all O2 to form FeO
                pSource[RealData::mFeO]   = mdotO2d/nO2FeO - mdotFeOevap;
                pSource[RealData::mFe3O4] = 0.0;
                epsilon = 0.70;
                dmO2dt = -mdotO2d;
                dmFeOformdt   = mdotO2d/nO2FeO;
                dmFe3O4formdt = 0;
            }
            pSource[RealData::mFe] = - pSource[RealData::mFeO]*nFeFeO - pSource[RealData::mFe3O4]*nFeFe3O4;
            pInt[IntData::regime] = 1;
            // std::cout << "diffusion-controlled" << std::endl;
        }
    }

    // Calculate Nusselt number, particle surface area, radiation parameters, gas conductivity, etc. for heat transfer
    Ap      = pi*dp*dp;
    SB      = 5.6704e-8;
    if (Tp > 1870){
        Nu = Nu*log(1.0+Bt)/Bt; // update Nusselt number with Stefan effect to use in heat transfer calculation
        epsilon = 0.7;
    }
    else {
        epsilon = 0.88;
    }
    conv    = -2*pi*rp*Nu*kgas*(Tp-Tgas);
    rad     = -Ap*epsilon*SB*(pow(Tp,4.0)-pow(Tgas,4.0));
    evap    = -hFeG(Tp)*mdotFe/M_Fe - hFeOG(Tp)*mdotFeO/M_Fe;
    dHpdt   = conv + rad + evap - dmO2dt*hO2(Tp)/M_O2;
    pSource[RealData::Hp] = dHpdt; 

    // std::cout << "kgas, nuStef" << kgas << " " << Nu << std::endl;
    // std::cout << "conv+rad, evap" << conv+rad << " " << evap << std::endl;
    // std::cout << "mass based hO2, addition" << hO2(Tp)/M_O2 << " " << dmO2dt*hO2(Tp)/M_O2 << std::endl;

    if (Tp > 1650){
        qFeO = qFeOl;
    }
    else{
        qFeO = qFeOs;
    }

    dHgOxidt = (mdotFe/M_Fe)*(0.5*qFe2O3s*M_Fe2O3) + (mdotFeO/M_FeO)*(0.5*qFe2O3s*M_Fe2O3-qFeOg*M_FeO); 

    // std::cout << "mFe: " << mFe << ", mFeO: " << mFeO << ", mFe3O4: " << mFe3O4 << std::endl;
    // std::cout << "dupdt: " << dupdt << ", dvpdt: " << dvpdt << std::endl;
    // std::cout << "conv loss: " << conv << ", rad loss: " << rad << ", mO2gain: " << -dmO2dt*hO2(Tp)/M_O2 << std::endl;
    // std::cout << "energy source for gas: " << -dHpdt << std::endl;

    qSource[gasVar::rho]    = dmO2dt;
    qSource[gasVar::rhou]   = -(mFe+mFeO+mFe3O4)*dupdt;
    if (spacedim == 1){
        qSource[gasVar::rhov]   = 0;
    }
    else {
        qSource[gasVar::rhov]   = -(mFe+mFeO+mFe3O4)*dvpdt;
    }
    qSource[gasVar::E]      = -dHpdt + dHgOxidt;  // + qFeO*dmFeOformdt + qFe3O4s*dmFe3O4formdt
    qSource[gasVar::rhoYO2] = dmO2dt;
    qSource[gasVar::rhoYN2] = 0;

    // Must assess phase transition here. If the particle temperature reaches a melting point temperature
    // (FeO:1650, Fe:1809, Fe3O4:1870 K) for the first time, we assign the particle enthalpy at that instance
    // to LFe, LFeO, LFe3O4, depending on which species is melting.
    if (Tp >= 1650){
        if (phaseFeO == 0){
            meltFeO = mFeO*HfuFeO;
            pSource[RealData::LFeO] += dHpdt + pSource[RealData::mFeO]*qFeO + pSource[RealData::mFe3O4]*qFe3O4s;
            if (fabs(pReal[RealData::LFeO])>meltFeO){
                pInt[IntData::FeO] = 1; 
                pReal[RealData::Hp] = Hparticle(mFe,mFeO,mFe3O4,1650.0,0,1,0);
            }
        }
        if (phaseFeO == 2){
            pInt[IntData::FeO] = 2;
            pReal[RealData::Hp] = Hparticle(mFe,mFeO,mFe3O4,1650.0,0,0,0);
        }
    }
    if (Tp >= 1809){
        if (phaseFe == 0){
            meltFe = mFe*HfuFe;
            pSource[RealData::LFe] += dHpdt + pSource[RealData::mFeO]*qFeO + pSource[RealData::mFe3O4]*qFe3O4s;
            if (fabs(pReal[RealData::LFe])>meltFe){
                pInt[IntData::Fe] = 1; 
                pReal[RealData::Hp] = Hparticle(mFe,mFeO,mFe3O4,1809,1,1,0);
            }
        }
        if (phaseFe == 2){
            pInt[IntData::Fe] = 2;
            pReal[RealData::Hp] = Hparticle(mFe,mFeO,mFe3O4,1809.0,0,1,0);
        }
    }
    if (Tp >= 1870){
        if (phaseFe3O4 == 0){
            meltFe3O4 = mFe3O4*HfuFe3O4;
            pSource[RealData::LFe3O4] += dHpdt + pSource[RealData::mFeO]*qFeO + pSource[RealData::mFe3O4]*qFe3O4s;
            if (fabs(pReal[RealData::LFe3O4])>meltFe3O4){
                pInt[IntData::Fe3O4] = 1; 
                pReal[RealData::Hp] = Hparticle(mFe,mFeO,mFe3O4,1870,1,1,1);
            }
        }        
        if (phaseFe3O4 == 2){
            pInt[IntData::Fe3O4] = 2;
            pReal[RealData::Hp] = Hparticle(mFe,mFeO,mFe3O4,1870.0,1,1,0);
        }
    }   
}

void
AmrLevelAdv::printParticleInfo()
{
  const int lev = 0;
  // const Geometry& geom = Geom(lev);
  // const Real* dx = geom.CellSize();
  // const Real* plo = geom.ProbLo();
  Real x, y=0, up, vp, mFe, mFeO, mFe3O4, Hp, Tp, LFe, LFeO, LFe3O4;
  int  phaseFe, phaseFeO, phaseFe3O4;

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){

    const int grid_id = mfi.index();
    const int tile_id = mfi.LocalTileIndex();
    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    auto& particles = particle_tile.GetArrayOfStructs();
    const int np = particles.numParticles();
    // std::cout << "number of particles in grid in printParticleInfo: " << np << std::endl;

    for(int pindex = 0; pindex < np; ++pindex) {
        ParticleType& p = particles[pindex];
        const IntVect& iv = this->Index(p, lev);
        x = p.pos(0);
        if (spacedim == 2){
            y = p.pos(1);
        }
        up     = p.rdata(RealData::u);
        vp     = p.rdata(RealData::v);
        mFe    = p.rdata(RealData::mFe);
        mFeO   = p.rdata(RealData::mFeO);
        mFe3O4 = p.rdata(RealData::mFe3O4);
        Hp     = p.rdata(RealData::Hp);
        LFe    = p.rdata(RealData::LFe);
        LFeO   = p.rdata(RealData::LFeO);
        LFe3O4 = p.rdata(RealData::LFe3O4);
        phaseFe    = p.idata(IntData::Fe);
        phaseFeO   = p.idata(IntData::FeO);
        phaseFe3O4 = p.idata(IntData::Fe3O4);

        Tp = Tparticle(mFe,mFeO,mFe3O4,Hp,phaseFe,phaseFeO,phaseFe3O4,LFe,LFeO,LFe3O4);
    }
  }
//   std::cout << "Position and velocity:: x: " << x << ", y: " << y << ", up: " << up << ", vp: " << vp << std::endl;
//   std::cout << "Masses:: Fe: " << mFe << ", FeO: " << mFeO << ", Fe3O4: " << mFe3O4 << ", mTot: " << mFe+mFeO+mFe3O4 << std::endl;
//   std::cout << "Enthalpy and Temperature:: Hp: " << Hp << ", Tp: " << Tp << std::endl;
}

void particleInit(double& energy0){

    int    phaseFe=0, phaseFeO=0, phaseFe3O4=0;
    // double deltaFeO, deltaFe3O4, rp0, rFeO0, rFe0;

    double Tp = TpInitial;

    // rp0        = 0.5*dp0;
    // deltaFeO   = 0.95*delta0;
    // deltaFe3O4 = 0.05*delta0;

    // rFeO0      = rp0*(1-deltaFe3O4);
    // rFe0       = rp0*(1-delta0);
    // mFe0       = rhoFe*(4.0/3.0)*pi*pow(rFe0,3.0);
    // mFeO0      = rhoFeO*(4.0/3.0)*pi*(pow(rFeO0,3.0)-pow(rFe0,3.0));
    // mFe3O40    = rhoFe3O4*(4.0/3.0)*pi*(pow(rp0,3.0)-pow(rFeO0,3.0));

    // interDist  = pow(1.0e3*(mFe0+mFeO0+mFe3O40)/1100.0,1.0/3.0);

    energy0    = Hparticle(mFe0,mFeO0,mFe3O40,Tp,phaseFe,phaseFeO,phaseFe3O4);

    // check mass values    
    // std::cout << "Fe, FeO, Fe3O4 mass: " << mFe0 << " " << mFeO0 << " " << mFe3O40 << ", energy: " << energy0 << std::endl;
}

Vector<double> getOxidationRates(const double& mFe, const double& mFeO, const double& mFe3O4, \
                                 const double& Tp,  const double& rp){

    Vector<double> rates(3);
    
    // if (Tp >= 1350){
    //     Abort("Particle ignited"); 
    // }

    Real vFe, vFeO, vFe3O4, rFe, rFeO, rFe3O4, XFeO, XFe3O4;
    Real FeOform, Fe3O4form, dmFedt, dmFeOdt, dmFe3O4dt;
    
    vFe = mFe/rhoFe;
    vFeO = mFeO/rhoFeO;
    vFe3O4 = mFe3O4/rhoFe3O4;

    rFe = pow(3.0*vFe/(4.0*pi),1.0/3.0);
    rFeO = pow(3.0*vFeO/(4.0*pi)+pow(rFe,3.0),1.0/3.0);
    rFe3O4 = pow(3.0*vFe3O4/(4.0*pi)+pow(rFeO,3.0),1.0/3.0);

    XFeO = rFeO-rFe;
    XFe3O4 = rFe3O4-rFeO;

    // std::cout << "XFeO, XFe3O4: " << XFeO << " " << XFe3O4 << std::endl; 

    FeOform = 4.0*pi*rhoFeO*rFeO*rFe*(1.0/XFeO)*k0FeOs*exp(-TaFeOs/Tp);
    Fe3O4form = 4.0*pi*rhoFe3O4*rFe3O4*rFeO*(1.0/XFe3O4)*k0Fe3O4s*exp(-TaFe3O4s/Tp);
    
    dmFedt = -FeOform*nFeFeO - Fe3O4form*nFeFe3O4;
    dmFeOdt = FeOform;
    dmFe3O4dt = Fe3O4form;

    rates[0] = dmFedt;
    rates[1] = dmFeOdt;
    rates[2] = dmFe3O4dt;

    return rates;
}

double Cdrag(double Re){
    if (Re < 1e-12){
        Re = 1e-12;
    }
    double Cd = (24.0/Re)*(1.0+0.15*pow(Re,0.678));
    return Cd;
}

double Hparticle(const double& mFe, const double& mFeO, const double& mFe3O4, const double& Tp, \
                 const int& phaseFe, const int& phaseFeO, const int& phaseFe3O4){

    // This function takes in the masses of each phase, flags of whether they are in solid- or liquid-phase,
    // and the particle temperature, in order to calculate the total enthalpy of the particle. 
    // The total enthalpy is returned in units of J.

    double eFe,eFeO,eFe3O4,eTotal;

    // if (Tp < 1809){
    //     eFe = hFeS(Tp);
    // }
    // else{
        if ((phaseFe == 0)||(phaseFe == 2)){ // Fe solid 
            eFe = hFeS(Tp);
        }
        else{ // Fe liquid
            eFe = hFeL(Tp);
        }
    // }
    // if (Tp < 1650){
    //     eFeO = hFeOS(Tp);
    // }
    // else{ 
        if ((phaseFeO == 0)||(phaseFeO == 2)){ // FeO solid
            eFeO = hFeOS(Tp);
        }
        else{ // FeO liquid
            eFeO = hFeOL(Tp);
        }
    // }
    // if (Tp < 1870){
    //     eFe3O4 = hFe3O4S(Tp);
    // }
    // else{ 
        if ((phaseFe3O4 == 0)||(phaseFe3O4 == 2)){ // Fe3O4 solid
            eFe3O4 = hFe3O4S(Tp);
        }
        else{ // Fe3O4 liquid
            eFe3O4 = hFe3O4L(Tp);
        }
    // }
    eTotal = mFe*(1/M_Fe)*eFe + mFeO*(1/M_FeO)*eFeO + mFe3O4*(1/M_Fe3O4)*eFe3O4;
    return eTotal;

}

// ____ TOTAL ENTHALPY POLYNOMIAL FITS ____ //


double hFeS(const double& Tp){
    double a1,a2,a3,a4,a5,a6,a7,b1,b2,hFeval;

    if (Tp < 500){
        a1 = 1.350490931E+04;
        a2 = -7.803806250E+02;
        a3 = 9.440171470E+00;
        a4 = -2.521767704E-02;
        a5 = 5.350170510E-05;
        a6 = -5.099094730E-08;
        a7 = 1.993862728E-11;
        b1 = 2.416521408E+03;
        b2 = -4.749002850E+01;
    }
    else if ((Tp >= 500) && (Tp < 800)){
       
        a1 = 3.543032740E+06;
        a2 = -2.447150531E+04;
        a3 = 6.561020930E+01;
        a4 = -7.043929680E-02;
        a5 = 3.181052870E-05;
        a6 = 0.000000000E+00;
        a7 = 0.000000000E+00;
        b1 = 1.345059978E+05;
        b2 = -4.133788690E+02;


    }
    else if ((Tp >= 800) && (Tp < 1042)){
        								
        a1 = 2.661026334E+09;
        a2 = -7.846827970E+06;
        a3 = -7.289212280E+02;
        a4 = 2.613888297E+01;
        a5 = -3.494742140E-02;
        a6 = 1.763752622E-05;
        a7 = -2.907723254E-09;
        b1 = 5.234868470E+07;
        b2 = -1.529052200E+04;
    }
    else if ((Tp >= 1042) && (Tp < 1184)){
        								
        a1 = 2.481923052E+08;
        a2 = 0.000000000E+00;
        a3 = -5.594349090E+02;
        a4 = 3.271704940E-01;
        a5 = 0.000000000E+00;
        a6 = 0.000000000E+00;
        a7 = 0.000000000E+00;
        b1 = 6.467503430E+05;
        b2 = 3.669168720E+03;
    }
    else if ((Tp >= 1184) && (Tp < 1665)){
        								
        a1 = 1.442428576E+09;
        a2 = -5.335491340E+06;
        a3 = 8.052828000E+03-0.091421931582384;
        a4 = -6.303089630E+00;
        a5 = 2.677273007E-03;
        a6 = -5.750045530E-07;
        a7 = 4.718611960E-11;
        b1 = 3.264264250E+07;
        b2 = -5.508852170E+04;
    }
    else {
        					
        a1 = -3.450190030E+08;
        a2 = 0.000000000E+00;
        a3 = 7.057501520E+02-0.151885932448342;
        a4 = -5.442977890E-01;
        a5 = 1.190040139E-04;
        a6 = 0.000000000E+00;
        a7 = 0.000000000E+00;
        b1 = -8.045725750E+05;
        b2 = -4.545180320E+03;
    }

    hFeval = R*Tp*(-a1/(Tp*Tp) + a2*log(Tp)/Tp + a3 + a4*Tp/2.0 + a5*Tp*Tp/3.0 + a6*Tp*Tp*Tp/4.0 + a7*Tp*Tp*Tp*Tp/5.0 + b1/Tp);
    return hFeval;    

}

double hFeL(const double& Tp){

    // NASA polynomial fits for total enthalpy of Fe in liquid-phase
    // Returns hFe in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,hFeval;
    a1 = 0.000000000e+00; 
    a2 = 0.000000000e+00;
    a3 = 5.535383324e+00; 
    a4 = 0.000000000e+00;
    a5 = 0.000000000e+00;
    a6 = 0.000000000e+00; 
    a7 = 0.000000000e+00; 
    b1 = -1.270608703e+03;
    b2 = -2.948115042e+01;

    hFeval = R*Tp*(-a1/(Tp*Tp) + a2*log(Tp)/Tp + a3 + a4*Tp/2.0 + a5*Tp*Tp/3.0 + a6*Tp*Tp*Tp/4.0 + a7*Tp*Tp*Tp*Tp/5.0 + b1/Tp);
    return hFeval;    

}

double hFeG(const double& Tp){

    // NASA polynomial fits for total enthalpy of Fe in gas-phase
    // Returns hFe in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,hFeval;
    if (Tp < 1000){
        a1=	6.790822660E+04;
        a2=	-1.197218407E+03;
        a3=	9.843393310E+00;
        a4=	-1.652324828E-02;
        a5=	1.917939959E-05;
        a6=	-1.149825371E-08;
        a7=	2.832773807E-12;
        b1=	5.466995940E+04;
        b2=	-3.383946260E+01;


    }
    else if ((Tp >= 1000) && (Tp < 6000)){
       
        a1=	-1.954923682E+06;
        a2=	6.737161100E+03;
        a3=	-5.486410970E+00;
        a4=	4.378803450E-03;
        a5=	-1.116286672E-06;
        a6=	1.544348856E-10;
        a7=	-8.023578182E-15;
        b1=	7.137370060E+03;
        b2=	6.504979860E+01;

    }
    else{
        								
        a1=	1.216352511E+09;
        a2=	-5.828563930E+05;
        a3=	9.789634510E+01;
        a4=	-5.370704430E-03;
        a5=	3.192037920E-08;
        a6=	6.267671430E-12;
        a7=	-1.480574914E-16;
        b1=	4.847648290E+06;
        b2=	-8.697289770E+02;

    }

    hFeval = R*Tp*(-a1/(Tp*Tp) + a2*log(Tp)/Tp + a3 + a4*Tp/2.0 + a5*Tp*Tp/3.0 + a6*Tp*Tp*Tp/4.0 + a7*Tp*Tp*Tp*Tp/5.0 + b1/Tp);
    return hFeval;    

}


double hFeOS(const double& Tp){

    // NASA polynomial fits for total enthalpy of FeO in solid-phase, for temperature-dependent phases
    // Returns hFeO in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,hFeOval;

    a1 = -1.179193966e+04; 
    a2 = 1.388393372e+02; 
    a3 = 2.999841854e+00; 
    a4 = 1.274527210e-02;
    a5 = -1.883886065e-05;
    a6 = 1.274258345e-08;
    a7 = -3.042206479e-12;
    b1 = -3.417350500e+04;
    b2 = -1.284759120e+01;

    hFeOval = R*Tp*(-a1/(Tp*Tp) + a2*log(Tp)/Tp + a3 + a4*Tp/2.0 + a5*Tp*Tp/3.0 + a6*Tp*Tp*Tp/4.0 + a7*Tp*Tp*Tp*Tp/5.0 + b1/Tp);
    return hFeOval;    

}

double hFeOL(const double& Tp){

    // NASA polynomial fits for total enthalpy of FeO in liquid-phase
    // Returns hFeO in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,hFeOval;
    a1 = 0.000000000e+00; 
    a2 = 0.000000000e+00;
    a3 = 8.147077819e+00; 
    a4 = 0.000000000e+00;
    a5 = 0.000000000e+00;
    a6 = 0.000000000e+00; 
    a7 = 0.000000000e+00; 
    b1 = -3.255080650e+04;
    b2 = -3.995344357e+01;
    
    hFeOval = R*Tp*(-a1/(Tp*Tp) + a2*log(Tp)/Tp + a3 + a4*Tp/2.0 + a5*Tp*Tp/3.0 + a6*Tp*Tp*Tp/4.0 + a7*Tp*Tp*Tp*Tp/5.0 + b1/Tp);
    return hFeOval;    

}

double hFeOG(const double& Tp){

    // NASA polynomial fits for total enthalpy of FeO in gas-phase
    // Returns hFe in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,hFeval;
    if (Tp < 1000){
        a1=	1.569282213E+04;
        a2=	-6.460188880E+01;
        a3=	2.458925470E+00;
        a4=	7.016047360E-03;
        a5=	-1.021405947E-05;
        a6=	7.179297870E-09;
        a7=	-1.978966365E-12;
        b1=	2.964572665E+04;
        b2=	1.326115545E+01;

    }
    else{
       
        a1=	-1.195971480E+05;
        a2=	-3.624864780E+02;
        a3=	5.518880750E+00;
        a4=	-9.978856890E-04;
        a5=	4.376913830E-07;
        a6=	-6.790629460E-11;
        a7=	3.639292680E-15;
        b1=	3.037985806E+04;
        b2=	-3.633655420E+00;


    }
    hFeval = R*Tp*(-a1/(Tp*Tp) + a2*log(Tp)/Tp + a3 + a4*Tp/2.0 + a5*Tp*Tp/3.0 + a6*Tp*Tp*Tp/4.0 + a7*Tp*Tp*Tp*Tp/5.0 + b1/Tp);
    return hFeval;    

}

double hFe3O4S(const double& Tp){

    // NASA polynomial fits for total enthalpy of Fe3O4 in solid-phase, for temperature-dependent phases
    // Returns hFe3O4 in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,hFe3O4val;

    if (Tp < 298){
        a1 = -5.182671230E+07;
        a2 = 1.293463453E+06;
        a3 = -1.341121962E+04;
        a4 = 7.401842440E+01;
        a5 = -2.285725885E-01;
        a6 = 3.752884300E-04;
        a7 = -2.559338368E-07;
        b1 = -5.570751240E+06;
        b2 = 6.575689710E+04;
    }
    else if ((Tp >= 298) && (Tp < 800)){

        a1 = -4.407671380E+06;
        a2 = 5.351760270E+04;
        a3 = -2.613667759E+02;
        a4 = 7.431931490E-01;
        a5 = -9.767843990E-04;
        a6 = 5.858865440E-07;
        a7 = -8.780843180E-11;
        b1 = -4.018075450E+05;
        b2 = 1.478276107E+03;

    }
    else if ((Tp >= 800) && (Tp < 850)){

        a1 = 0.000000000E+00;
        a2 = 0.000000000E+00;
        a3 = -1.070116148E+02;
        a4 = 1.738436706E-01;
        a5 = 0.000000000E+00;
        a6 = 0.000000000E+00;
        a7 = 0.000000000E+00;
        b1 = -9.231022600E+04;
        b2 = 6.169265720E+02;

    }
    else {

        a1 = 5.731691980e+07;
        a2 = -1.816105186e+05;
        a3 = 2.777813396e+02;
        a4 = -1.849830315e-01;
        a5 = 6.145641150e-05;
        a6 = -3.660350860e-09;
        a7 = -1.383713617e-12;
        b1 = 9.906992850e+05;
        b2 = -1.868855147E+03;
    }

    hFe3O4val = R*Tp*(-a1/(Tp*Tp) + a2*log(Tp)/Tp + a3 + a4*Tp/2.0 + a5*Tp*Tp/3.0 + a6*Tp*Tp*Tp/4.0 + a7*Tp*Tp*Tp*Tp/5.0 + b1/Tp);
    return hFe3O4val;    

}

double hFe3O4L(const double& Tp){

    // NASA polynomial fits for total enthalpy of Fe3O4 in liquid-phase
    // Returns hFeO in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,hFe3O4val;
    a1 = 0.000000000e+00; 
    a2 = 0.000000000e+00;
    a3 = 2.415439996e+01; 
    a4 = 0.000000000e+00;
    a5 = 0.000000000e+00;
    a6 = 0.000000000e+00; 
    a7 = 0.000000000e+00; 
    b1 = -1.242541349e+05;
    b2 = -1.109781572e+02;
    
    hFe3O4val = R*Tp*(-a1/(Tp*Tp) + a2*log(Tp)/Tp + a3 + a4*Tp/2.0 + a5*Tp*Tp/3.0 + a6*Tp*Tp*Tp/4.0 + a7*Tp*Tp*Tp*Tp/5.0 + b1/Tp);
    return hFe3O4val;    

}

// ____ SPECIFIC HEAT CAPACITY POLYNOMIAL FITS ____ //

double cpFeS(const double& Tp){

    // NASA polynomial fits for specific heat capacity of Fe in liquid-phase
    // Returns cpFe in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,cpFeval;

    if (Tp < 500){
        a1=	1.350490931E+04;
        a2=	-7.803806250E+02;
        a3=	9.440171470E+00;
        a4=	-2.521767704E-02;
        a5=	5.350170510E-05;
        a6=	-5.099094730E-08;
        a7=	1.993862728E-11;
        b1=	2.416521408E+03;
        b2=	-4.749002850E+01;

    }
    else if ((Tp >= 500) && (Tp < 800)){
       
        a1=	3.543032740E+06;
        a2=	-2.447150531E+04;
        a3=	6.561020930E+01;
        a4=	-7.043929680E-02;
        a5=	3.181052870E-05;
        a6=	0.000000000E+00;
        a7=	0.000000000E+00;
        b1=	1.345059978E+05;
        b2=	-4.133788690E+02;



    }
    else if ((Tp >= 800) && (Tp < 1042)){
        								
        a1=	2.661026334E+09;
        a2=	-7.846827970E+06;
        a3=	-7.289212280E+02;
        a4=	2.613888297E+01;
        a5=	-3.494742140E-02;
        a6=	1.763752622E-05;
        a7=	-2.907723254E-09;
        b1=	5.234868470E+07;
        b2=	-1.529052200E+04;

    }
    else if ((Tp >= 1042) && (Tp < 1184)){
        								
        a1=	2.481923052E+08;
        a2=	0.000000000E+00;
        a3=	-5.594349090E+02;
        a4=	3.271704940E-01;
        a5=	0.000000000E+00;
        a6=	0.000000000E+00;
        a7=	0.000000000E+00;
        b1=	6.467503430E+05;
        b2=	3.669168720E+03;

    }
    else if ((Tp >= 1184) && (Tp < 1665)){
        								
        a1=	 1.442428576E+09;
        a2=	-5.335491340E+06;
        a3=	 8.052828000E+03+0.905406354181551;
        a4=	-6.303089630E+00;
        a5=	 2.677273007E-03;
        a6=	-5.750045530E-07;
        a7=	 4.718611960E-11;
        b1=	 3.264264250E+07;
        b2=	-5.508852170E+04;

    }
    else {
        					
        a1=	-3.450190030E+08;
        a2=	0.000000000E+00;
        a3=	7.057501520E+02+0.519144362184648;
        a4=	-5.442977890E-01;
        a5=	1.190040139E-04;
        a6=	0.000000000E+00;
        a7=	0.000000000E+00;
        b1=	-8.045725750E+05;
        b2=	-4.545180320E+03;

    }

    cpFeval = R*(a1/(Tp*Tp) + a2/Tp + a3 + a4*Tp + a5*Tp*Tp + a6*Tp*Tp*Tp + a7*Tp*Tp*Tp*Tp);
    return cpFeval;    

}

double cpFeL(const double& Tp){

    // NASA polynomial fits for specific heat capacity of Fe in liquid-phase
    // Returns cpFe in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,cpFeval;
    a1 = 0.000000000e+00; 
    a2 = 0.000000000e+00;
    a3 = 5.535383324e+00; 
    a4 = 0.000000000e+00;
    a5 = 0.000000000e+00;
    a6 = 0.000000000e+00; 
    a7 = 0.000000000e+00; 
    b1 = -1.270608703e+03;
    b2 = -2.948115042e+01;


    cpFeval = R*(a1/(Tp*Tp) + a2/Tp + a3 + a4*Tp + a5*Tp*Tp + a6*Tp*Tp*Tp + a7*Tp*Tp*Tp*Tp);
    return cpFeval;    

}

double cpFeG(const double& Tp){

    // NASA polynomial fits for specific heat capacity of Fe in gas-phase
    // Returns cpFe in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,cpFeval;
    if (Tp < 1000){
        a1=	6.790822660E+04;
        a2=	-1.197218407E+03;
        a3=	9.843393310E+00;
        a4=	-1.652324828E-02;
        a5=	1.917939959E-05;
        a6=	-1.149825371E-08;
        a7=	2.832773807E-12;
        b1=	5.466995940E+04;
        b2=	-3.383946260E+01;
    }
    else if ((Tp >= 1000) && (Tp < 6000)){
       
        a1=	-1.954923682E+06;
        a2=	6.737161100E+03;
        a3=	-5.486410970E+00;
        a4=	4.378803450E-03;
        a5=	-1.116286672E-06;
        a6=	1.544348856E-10;
        a7=	-8.023578182E-15;
        b1=	7.137370060E+03;
        b2=	6.504979860E+01;

    }
    else{
        								
        a1=	1.216352511E+09;
        a2=	-5.828563930E+05;
        a3=	9.789634510E+01;
        a4=	-5.370704430E-03;
        a5=	3.192037920E-08;
        a6=	6.267671430E-12;
        a7=	-1.480574914E-16;
        b1=	4.847648290E+06;
        b2=	-8.697289770E+02;

    }

    cpFeval = R*(a1/(Tp*Tp) + a2/Tp + a3 + a4*Tp + a5*Tp*Tp + a6*Tp*Tp*Tp + a7*Tp*Tp*Tp*Tp);
    return cpFeval;

}

double cpFeOS(const double& Tp){

    // NASA polynomial fits for specific heat capacity of FeO in solid-phase, for temperature-dependent phases
    // Returns cpFeO in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,cpFeOval;

    a1=	-1.179193966e+04;
    a2=	1.388393372e+02;
    a3=	2.999841854e+00;
    a4=	1.274527210e-02;
    a5=	-1.883886065e-05;
    a6=	1.274258345e-08;
    a7=	-3.042206479e-12;
    b1=	-3.417350500e+04;
    b2=	-1.284759120E+01;


    cpFeOval = R*(a1/(Tp*Tp) + a2/Tp + a3 + a4*Tp + a5*Tp*Tp + a6*Tp*Tp*Tp + a7*Tp*Tp*Tp*Tp);
    return cpFeOval;    

}

double cpFeOL(const double& Tp){

    // NASA polynomial fits for specific heat capacity of FeO in liquid-phase
    // Returns cpFeO in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,cpFeOval;
    a1 = 0.000000000e+00; 
    a2 = 0.000000000e+00;
    a3 = 8.147077819e+00; 
    a4 = 0.000000000e+00;
    a5 = 0.000000000e+00;
    a6 = 0.000000000e+00; 
    a7 = 0.000000000e+00; 
    b1 = -3.255080650e+04;
    b2 = -3.995344357e+01;
    
    cpFeOval = R*(a1/(Tp*Tp) + a2/Tp + a3 + a4*Tp + a5*Tp*Tp + a6*Tp*Tp*Tp + a7*Tp*Tp*Tp*Tp);
    return cpFeOval;    

}

double cpFeOG(const double& Tp){

    // NASA polynomial fits for specific heat capacity of FeO in gas-phase
    // Returns cpFeO in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,cpFeOval;
    if (Tp < 1000){
        a1=	1.569282213E+04;
        a2=	-6.460188880E+01;
        a3=	2.458925470E+00;
        a4=	7.016047360E-03;
        a5=	-1.021405947E-05;
        a6=	7.179297870E-09;
        a7=	-1.978966365E-12;
        b1=	2.964572665E+04;
        b2=	1.326115545E+01;
    }
    else{
       
        a1=	-1.195971480E+05;
        a2=	-3.624864780E+02;
        a3=	5.518880750E+00;
        a4=	-9.978856890E-04;
        a5=	4.376913830E-07;
        a6=	-6.790629460E-11;
        a7=	3.639292680E-15;
        b1=	3.037985806E+04;
        b2=	-3.633655420E+00;
    }

    cpFeOval = R*(a1/(Tp*Tp) + a2/Tp + a3 + a4*Tp + a5*Tp*Tp + a6*Tp*Tp*Tp + a7*Tp*Tp*Tp*Tp);
    return cpFeOval;

}

double cpFe3O4S(const double& Tp){

    // NASA polynomial fits for specific heat capacity of Fe3O4 in solid-phase, for temperature-dependent phases
    // Returns cpFe3O4 in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,cpFe3O4val;

    if (Tp < 298){
        a1=	-5.182671230E+07;
        a2=	1.293463453E+06;
        a3=	-1.341121962E+04;
        a4=	7.401842440E+01;
        a5=	-2.285725885E-01;
        a6=	3.752884300E-04;
        a7=	-2.559338368E-07;
        b1=	-5.570751240E+06;
        b2=	6.575689710E+04;

    }
    else if ((Tp >= 298) && (Tp < 800)){

        a1=	-4.407671380E+06;
        a2=	5.351760270E+04;
        a3=	-2.613667759E+02;
        a4=	7.431931490E-01;
        a5=	-9.767843990E-04;
        a6=	5.858865440E-07;
        a7=	-8.780843180E-11;
        b1=	-4.018075450E+05;
        b2=	1.478276107E+03;


    }
    else if ((Tp >= 800) && (Tp < 850)){

        a1=	0.000000000E+00;
        a2=	0.000000000E+00;
        a3=	-1.070116148E+02;
        a4=	1.738436706E-01;
        a5=	0.000000000E+00;
        a6=	0.000000000E+00;
        a7=	0.000000000E+00;
        b1=	-9.231022600E+04;
        b2=	6.169265720E+02;


    }
    else {

        a1=	5.731691980e+07;
        a2=	-1.816105186e+05;
        a3=	2.777813396e+02;
        a4=	-1.849830315e-01;
        a5=	6.145641150e-05;
        a6=	-3.660350860e-09;
        a7=	-1.383713617e-12;
        b1=	9.906992850e+05;
        b2=	-1.868855147E+03;

    }

    cpFe3O4val = R*(a1/(Tp*Tp) + a2/Tp + a3 + a4*Tp + a5*Tp*Tp + a6*Tp*Tp*Tp + a7*Tp*Tp*Tp*Tp);
    return cpFe3O4val;    

}

double cpFe3O4L(const double& Tp){

    // NASA polynomial fits for specific heat capacity of Fe3O4 in liquid-phase
    // Returns cpFeO in J/mol

    double a1,a2,a3,a4,a5,a6,a7,b1,b2,cpFe3O4val;
    a1 = 0.000000000e+00; 
    a2 = 0.000000000e+00;
    a3 = 2.415439996e+01; 
    a4 = 0.000000000e+00;
    a5 = 0.000000000e+00;
    a6 = 0.000000000e+00; 
    a7 = 0.000000000e+00; 
    b1 = -1.242541349e+05;
    b2 = -1.109781572e+02;
    
    cpFe3O4val = R*(a1/(Tp*Tp) + a2/Tp + a3 + a4*Tp + a5*Tp*Tp + a6*Tp*Tp*Tp + a7*Tp*Tp*Tp*Tp);
    return cpFe3O4val;    

}

double cpparticle(const double& mFe, const double& mFeO, const double& mFe3O4, const double& Tp, \
                  const int& phaseFe, const int& phaseFeO, const int& phaseFe3O4){

    // This function takes in the masses of each phase, flags of whether they are in solid- or liquid-phase,
    // and the particle temperature, in order to calculate the total heat capacity of the particle. 
    // The heat capacity is returned in units of J/K.

    double cpFe,cpFeO,cpFe3O4,cpTotal;

    // if (Tp < 1809){
    //     cpFe = cpFeS(Tp);
    // }
    // else{
        if ((phaseFe == 0)||(phaseFe == 2)){ // Fe solid 
            cpFe = cpFeS(Tp);
        }
        else{ // Fe liquid
            cpFe = cpFeL(Tp);
        }
    // }
    // if (Tp < 1650){
    //     cpFeO = cpFeOS(Tp);
    // }
    // else{
        if ((phaseFeO == 0)||(phaseFeO == 2)) { // FeO solid
            cpFeO = cpFeOS(Tp);
        }
        else{ // FeO liquid
            cpFeO = cpFeOL(Tp);
        }
    // }
    // if (Tp < 1870){
    //     cpFe3O4 = cpFe3O4S(Tp);
    // }
    // else{
        if ((phaseFe3O4 == 0)||(phaseFe3O4 == 2)){ // Fe3O4 solid
            cpFe3O4 = cpFe3O4S(Tp);
        }
        else{ // Fe3O4 liquid
            cpFe3O4 = cpFe3O4L(Tp);
        }
    // }

    cpTotal = mFe*(1/M_Fe)*cpFe+mFeO*(1/M_FeO)*cpFeO+mFe3O4*(1/M_Fe3O4)*cpFe3O4;
    return cpTotal;

}

double Tparticle(const double& mFe, const double& mFeO, const double& mFe3O4, const double& Hp, \
                 int& phaseFe, int& phaseFeO, int& phaseFe3O4, \
                 const double& LFe, const double& LFeO, const double& LFe3O4){

    // This function takes in the phase masses, solid- or liquid-phase flags, particle temperature,
    // and the phase transition progress variables (LFe, LFeO, LFe3O4) to iteratively determine
    // the particle temperature.

    // If a phase flag is still solid but the progress variable is non-zero, we set the particle
    // temperature to the phase transition temperature of that phase.

    double Tpn;

    if ((phaseFe == 0)&&(LFe != 0.0)){
        Tpn = 1809;
        // std::cout << "in Fe melt, Tp = 1809" << std::endl;
        return Tpn;
    }
    if ((phaseFeO == 0)&&(LFeO != 0.0)){
        Tpn = 1650;
        // std::cout << "in FeO melt, Tp = 1650" << std::endl;
        return Tpn;
    }
    if ((phaseFe3O4 == 0)&&(LFe3O4 != 0.0)){
        Tpn = 1870;
        // std::cout << "in Fe3O4 melt, Tp = 1870" << std::endl;
        return Tpn;
    }

    // Otherwise, we shall proceed with the familiar fixed point algorithm to solve for the particle temperature
    // based on the true enthalpy value, and the temperature-dependent polynomial fits for enthalpy and heat capacity
    // of each species in the particle.

    double Tp0=300,tol=1e-4,error=1;
    double dH,cp,iter=0;

    // std::cout << "particle enthalpy " << Hp << std::endl;

    while (error > tol)
    {
        dH    = Hp - Hparticle(mFe, mFeO, mFe3O4, Tp0, phaseFe, phaseFeO, phaseFe3O4);
        cp    = cpparticle(mFe, mFeO, mFe3O4, Tp0, phaseFe, phaseFeO, phaseFe3O4);
        Tpn   = Tp0 + dH/cp;
        error = std::fabs(Tpn-Tp0)/std::fabs(Tpn);
        Tp0 = Tpn;
        iter += 1;
        // std::cout << "Tpn: " << Tpn << std::endl;
        if (iter > 20){
            Abort("Particle temperature iteration does not converge"); 
        }
    }
    // Conditional statements for solid-to-liquid phase transition
    if ((Tpn >= 1650)&&(phaseFeO == 0)){
        // std::cout << "in FeO melt, Tp = 1650" << std::endl;
        Tpn = 1650;
        return Tpn;
    }
    if ((Tpn >= 1809)&&(phaseFe == 0)){
        // std::cout << "in Fe melt, Tp = 1809" << std::endl;
        Tpn = 1809;
        return Tpn;
    }
    if ((Tpn >= 1870)&&(phaseFe3O4 == 0)){
        // std::cout << "in Fe3O4 melt, Tp = 1870" << std::endl;
        Tpn = 1870;
        return Tpn;
    }
    // Conditional statements for solid-to-liquid phase transition
    if ((Tpn <= 1650)&&(phaseFeO == 1)){
        // std::cout << "in FeO melt, Tp = 1650" << std::endl;
        Tpn = 1650;
        phaseFeO = 2;
        return Tpn;
    }
    if ((Tpn <= 1809)&&(phaseFe == 1)){
        // std::cout << "in Fe melt, Tp = 1809" << std::endl;
        Tpn = 1809;
        phaseFe = 2;
        return Tpn;
    }
    if ((Tpn <= 1870)&&(phaseFe3O4 == 1)){
        // std::cout << "in Fe3O4 melt, Tp = 1870" << std::endl;
        Tpn = 1870;
        phaseFe3O4 = 2;
        return Tpn;
    }
    // std::cout << "Particle temperature: " << Tpn << std::endl;
    return Tpn;

}

double filmAverage(const double& Tp, const double& Tg){
    double fraction = 1.0/2.0;
    double Tfilm = fraction*Tg + (1.0-fraction)*Tp;
    return Tfilm;
}
