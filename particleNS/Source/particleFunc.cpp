#include "eulerFunc.H"
#include "diffusionFunc.H"
#include "thermoTransport.H"
#include "particleFunc.H"
#include "AmrLevelAdv.H"
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
extern const int spacedim;
extern double R;
extern double pi;
extern double dp0;               // mu m                 Initial particle diameter
extern double delta0;            // [-]                  Initial oxide layer thickness to particle size ratio
extern double rhoFe;             // kg/m^3               Solid-phase iron(Fe) density
extern double rhoFeO;            // kg/m^3               Solid-phase FeO density
extern double rhoFe3O4;          // kg/m^3               Solid-phase Fe3O4 density
extern double M_Fe;              // kg/mol               Molar mass of iron (Fe)
extern double M_FeO;             // kg/mol               Molar mass of wustite (FeO)
extern double M_Fe3O4;           // kg/mol               Molar mass of wustite (FeO)
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
extern double M_O2;

using namespace amrex;

// Use this file to write functions required for particle calculations

void
AmrLevelAdv::initParticles ()
{
  const int lev = 0;
  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
  {
      auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                      mfi.LocalTileIndex())];

      ParticleType p;
      p.id()   = ParticleType::NextID();
      p.cpu()  = ParallelDescriptor::MyProc();

      p.pos(0) = 0.0;
      // p.pos(1) =  0.0;

      double mFe0,mFeO0,mFe3O40,energy0;
      particleInit(mFe0,mFeO0,mFe3O40,energy0);

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

      particles.push_back(p);
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

  std::cout << "in getParticleInfo" << std::endl;

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){

    const int grid_id = mfi.index();
    const int tile_id = mfi.LocalTileIndex();
    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    auto& particles = particle_tile.GetArrayOfStructs();
    const int np = particles.numParticles();
    std::cout << "number of particles in grid in getParticleInfo, within advance: " << np << std::endl;

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
  Real x, y, lx, ly;
  int  i=0, j=0, k=0;
  Vector<double> q(NUM_STATE);
  Vector<double> qSource(NUM_STATE,0);
  Vector<double> pReal(RealData::ncomps);
  Vector<int>    pInt(IntData::ncomps);
  Vector<double> pSource(RealData::ncomps,0);

  std::cout << "in updateParticleInfo" << std::endl;

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){

    const int grid_id = mfi.index();
    const int tile_id = mfi.LocalTileIndex();
    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    auto& particles = particle_tile.GetArrayOfStructs();
    const int np = particles.numParticles();
    std::cout << "number of particles in grid in getParticleInfo, within advance: " << np << std::endl;

    // Access S_new multifab information here
    const Box& bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);
    const auto& arr = Sborder.array(mfi);

    for(int pindex = 0; pindex < np; ++pindex) {
        ParticleType& p = particles[pindex];
        const IntVect& iv = this->Index(p, lev);

        x  = p.pos(0);
        lx = (x - lo.x) / dx;
        i  = static_cast<int>(Math::floor(lx));
        if (spacedim == 2){
            y  = p.pos(1);
            ly = (y - lo.y) / dy;
            j  = static_cast<int>(Math::floor(ly));
        }

        std::cout << "position in i is: " << lx << ", cell number is: " << i << std::endl;
        std::cout << "density in i is: " << arr(i,j,k,0) << std::endl;

        for (int h = 0; h < NUM_STATE; h++){
            q[h] = arr(i,j,k,h);
        }
        for (int h = 0; h < RealData::ncomps; h++){
            pReal[h] = p.rdata(h);
        }
        for (int h = 0; h < IntData::ncomps; h++){
            pInt[h] = p.idata(h);
        }

        getSource(qSource,pSource,arr,i,j,k,pReal,pInt,dt,dx,dy);

        p.pos(0) = p.pos(0) + dt*p.rdata(RealData::u);
        if (spacedim == 2){
            p.pos(1) = p.pos(1) + dt*p.rdata(RealData::v);
        }
        for (int h = 0; h < RealData::ncomps; h++){
            p.rdata(h) = p.rdata(h) + dt*pSource[h];
        }
    }
  }
  Redistribute();
}

void getSource(Vector<double>& qSource, Vector<double>& pSource, auto& arr, \
               const int& i, const int& j, const int& k, const Vector<double>& pReal, \
               const Vector<int>& pInt, const double& dt, const double& dx, const double& dy){
    
    // Declare particle and gas variables
    Real up, vp, wp, mFe, mFeO, mFe3O4, Hp, LFe, LFeO, LFe3O4;
    int  phaseFe, phaseFeO, phaseFe3O4;
    Real rho, rhou, rhov, energy, rhoYO2, rhoYN2, u, v, p, Tgas, YO2, YN2;

    // Declare variables used for drag calculation
    Real Tp, vTot, rp, dp, rhop, Re, mu, uRel, vRel, CD;

    // Declare variables used for heat transfer calculations
    Real Nu, Ap, epsilon, SB, kgas, conv, rad;


    // Declare source terms
    Real dupdt, dvpdt, dHpdt, dmO2dt; //, dHpdt, dmFedt, dmFeOdt, dmFe3O4dt;
    Vector<double> dmdt(3);

    rho    = arr(i,j,k,0);
    rhou   = arr(i,j,k,1);
    rhov   = arr(i,j,k,2);
    energy = arr(i,j,k,3);
    rhoYO2 = arr(i,j,k,4);
    rhoYN2 = arr(i,j,k,5);

    u    = rhou/rho;
    v    = rhov/rho;
    YO2  = rhoYO2/rho;
    YN2  = rhoYN2/rho;
    Tgas = Tg(rho,u,v,YO2,YN2,energy);
    p    = pressure(rho,YO2,YN2,Tgas);

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

    // Calculate necessary particle parameters here:
    vTot = mFe/rhoFe+mFeO/rhoFeO+mFe3O4/rhoFe3O4;    // total particle volume, m^3
    rp   = pow(3.0*vTot*0.25/pi,1.0/3.0);            // particle outer radius, m
    dp   = 2*rp;
    rhop = (mFe+mFeO+mFe3O4)/vTot;                   // average particle density, kg/m^3
    
    // Calculate Reynolds number for particle flow here:
    uRel = fabs(u-up);
    mu = muMix(muO2(Tgas),muN2(Tgas),YO2,YN2);
    Re = 2.0*rp*rho*uRel/mu;

    // Calculate drag coefficient here:
    CD = Cdrag(Re);

    // Calculate the rate of change of particle velocity here:
    dupdt = (3*CD*rho/(8*rp*rhop))*(u-up)*uRel;
    pSource[RealData::u] = dupdt; 

    if (spacedim == 2){
        vRel  = fabs(v-vp);
        Re    = 2.0*rp*rho*vRel/mu;
        CD    = Cdrag(Re);
        dvpdt = (3*CD*rho/(8*rp*rhop))*(v-vp)*vRel;
        pSource[RealData::v] = dvpdt; 
    }

    // Calculate change of Fe, FeO, and Fe3O4 mass based on temperature and oxide layer
    dmdt = getOxidationRates(mFe,mFeO,mFe3O4,Tp,rp);
    pSource[RealData::mFe] = dmdt[0];
    pSource[RealData::mFeO] = dmdt[1];
    pSource[RealData::mFe3O4] = dmdt[2];

    // Get rate of oxygen mass consumption due to reactions
    dmO2dt = -( nO2FeO*dmdt[1] +  nO2Fe3O4*dmdt[2]);

    // Calculate Nusselt number, particle surface area, radiation parameters, gas conductivity, etc. for heat transfer
    Nu      = 2.0;
    Ap      = pi*dp*dp;
    epsilon = 0.88;
    SB      = 5.6704e-8;
    double Tfilm = (1.0/3.0)*(2*Tp+Tgas);
    kgas    = kMix(kO2(Tfilm), kN2(Tfilm), YO2, YN2);
    conv    = -Ap*kgas*Nu*(1.0/dp)*(Tp-Tgas);
    rad     = Ap*epsilon*SB*(pow(Tp,4.0)-pow(Tgas,4.0));
    dHpdt   = conv + rad - dmO2dt*hO2(Tp)/M_O2;
    pSource[RealData::Hp] = dHpdt; 

    
    // std::cout << "rp: " << rp << " uRel: " << uRel << " CD: " << CD << " rho: " << rho << " rhop: " << rhop << std::endl;

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
    std::cout << "number of particles in grid in printParticleInfo: " << np << std::endl;

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
  std::cout << "Position and velocity:: x: " << x << ", y: " << y << ", up: " << up << ", vp: " << vp << std::endl;
  std::cout << "Masses:: Fe: " << mFe << ", FeO: " << mFeO << ", Fe3O4: " << mFe3O4 << ", mTot: " << mFe+mFeO+mFe3O4 << std::endl;
  std::cout << "Enthalpy and Temperature:: Hp: " << Hp << ", Tp: " << Tp << std::endl;
}

void particleInit(double& mFe0, double& mFeO0, double& mFe3O40, double& energy0){

    int    phaseFe=0, phaseFeO=0, phaseFe3O4=0;
    double deltaFeO, deltaFe3O4, rp0, rFeO0, rFe0;

    double Tp = TpInitial;

    rp0        = 0.5*dp0;
    deltaFeO   = 0.95*delta0;
    deltaFe3O4 = 0.05*delta0;

    rFeO0      = rp0*(1-deltaFe3O4);
    rFe0       = rp0*(1-delta0);
    mFe0       = rhoFe*(4.0/3.0)*pi*pow(rFe0,3.0);
    mFeO0      = rhoFeO*(4.0/3.0)*pi*(pow(rFeO0,3.0)-pow(rFe0,3.0));
    mFe3O40    = rhoFe3O4*(4.0/3.0)*pi*(pow(rp0,3.0)-pow(rFeO0,3.0));

    energy0    = Hparticle(mFe0,mFeO0,mFe3O40,Tp,phaseFe,phaseFeO,phaseFe3O4);

    // check mass values
    std::cout << "Fe, FeO, Fe3O4 mass: " << mFe0 << " " << mFeO0 << " " << mFe3O40 << ", energy: " << energy0 << std::endl;
}

Vector<double> getOxidationRates(const double& mFe, const double& mFeO, const double& mFe3O4, \
                                 const double& Tp,  const double& rp){

    Vector<double> rates(3);
    
    if (Tp >= 1600){
        // for (int i = 0; i < 3; i++){
        //     rates[i] = 0;
        // }
        // return rates;
        Abort("Particle ignited"); 
    }

    Real vFe, vFeO, vFe3O4, rFe, rFeO, XFeO, XFe3O4;
    Real FeOform, Fe3O4form, dmFedt, dmFeOdt, dmFe3O4dt;
    
    vFe = mFe/rhoFe;
    vFeO = mFeO/rhoFeO;
    vFe3O4 = mFe3O4/rhoFe3O4;
    
    rFe = pow(3.0*vFe/(4.0*pi),1.0/3.0);
    rFeO = pow(3.0*vFeO/(4.0*pi)+pow(rFe,3.0),1.0/3.0);

    XFeO = rFeO-rFe;
    XFe3O4 = rp-rFeO;

    FeOform = 4*pi*rhoFeO*rFeO*rFe*(1/XFeO)*k0FeOs*exp(-TaFeOs/Tp);
    Fe3O4form = 4*pi*rhoFe3O4*rp*rFeO*(1/XFe3O4)*k0Fe3O4s*exp(-TaFe3O4s/Tp);

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

    if (phaseFe == 0){ // Fe solid 
        eFe = hFeS(Tp);
    }
    else{ // Fe liquid
        eFe = hFeL(Tp);
    }
    if (phaseFeO == 0){ // FeO solid
        eFeO = hFeOS(Tp);
    }
    else{ // FeO liquid
        eFeO = hFeOL(Tp);
    }
    if (phaseFe3O4 == 0){ // Fe3O4 solid
        eFe3O4 = hFe3O4S(Tp);
    }
    else{ // Fe3O4 liquid
        eFe3O4 = hFe3O4L(Tp);
    }

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
        a3 = -5.594349090E+02+0.0914;
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
        a3 = 8.052828000E+03;
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
        a3 = 7.057501520E+02;
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
        a3=	-5.594349090E+02+0.0914;
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
        a3=	 8.052828000E+03;
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
        a3=	7.057501520E+02;
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

    if (phaseFe == 0){ // Fe solid 
        cpFe = cpFeS(Tp);
    }
    else{ // Fe liquid
        cpFe = cpFeL(Tp);
    }
    if (phaseFeO == 0){ // FeO solid
        cpFeO = cpFeOS(Tp);
    }
    else{ // FeO liquid
        cpFeO = cpFeOL(Tp);
    }
    if (phaseFe3O4 == 0){ // Fe3O4 solid
        cpFe3O4 = cpFe3O4S(Tp);
    }
    else{ // Fe3O4 liquid
        cpFe3O4 = cpFe3O4L(Tp);
    }

    cpTotal = mFe*(1/M_Fe)*cpFe+mFeO*(1/M_FeO)*cpFeO+mFe3O4*(1/M_Fe3O4)*cpFe3O4;
    return cpTotal;

}

double Tparticle(const double& mFe, const double& mFeO, const double& mFe3O4, const double& Hp, \
                 const int& phaseFe, const int& phaseFeO, const int& phaseFe3O4, \
                 const double& LFe, const double& LFeO, const double& LFe3O4){

    // This function takes in the phase masses, solid- or liquid-phase flags, particle temperature,
    // and the phase transition progress variables (LFe, LFeO, LFe3O4) to iteratively determine
    // the particle temperature.

    // If a phase flag is still solid but the progress variable is non-zero, we set the particle
    // temperature to the phase transition temperature of that phase.

    double Tpn;

    if ((phaseFe == 0)&&(LFe != 0)){
        Tpn = 1809;
        return Tpn;
    }
    if ((phaseFeO == 0)&&(LFeO != 0)){
        Tpn = 1650;
        return Tpn;
    }
    if ((phaseFe3O4 == 0)&&(LFe3O4 != 0)){
        Tpn = 1870;
        return Tpn;
    }

    // Otherwise, we shall proceed with the familiar fixed point algorithm to solve for the particle temperature
    // based on the true enthalpy value, and the temperature-dependent polynomial fits for enthalpy and heat capacity
    // of each species in the particle.

    double Tp0=300,tol=1e-3,error=1;
    double dH,cp,iter=0;

    while (error > tol)
    {
        dH    = Hp - Hparticle(mFe, mFeO, mFe3O4, Tp0, phaseFe, phaseFeO, phaseFe3O4);
        cp    = cpparticle(mFe, mFeO, mFe3O4, Tp0, phaseFe, phaseFeO, phaseFe3O4);
        Tpn   = Tp0 + dH/cp;
        error = std::fabs(Tpn-Tp0)/Tpn;
        std::cout << "ener: " << Hp << ", enthalpy diff: " << dH << \
        ", Tpn: " << Tpn << ", error: " << error << std::endl;
        Tp0 = Tpn;
        iter += 1;
        if (iter > 20){
            Abort("Particle temperature iteration does not converge"); 
        }
    }

    return Tpn;

}
