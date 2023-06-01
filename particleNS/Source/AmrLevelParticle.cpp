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

using namespace amrex;

void
AmrLevelAdv::initParticles ()
{
  const int lev = 0;
  // using MyParIter = amrex::ParIter<RealData::ncomps, IntData::ncomps>;
  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
  {
      auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                      mfi.LocalTileIndex())];

      ParticleType p;
      p.id()   = ParticleType::NextID();
      p.cpu()  = ParallelDescriptor::MyProc();

      p.pos(0) = 0.5;
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

  std::cout << "particle has been initialized. to confirm, we search for the particle and get its info" << std::endl;

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){

      const int grid_id = mfi.index();
      const int tile_id = mfi.LocalTileIndex();
      auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
      auto& particles = particle_tile.GetArrayOfStructs();
      const int np = particles.numParticles();
      std::cout << "number of particles in grid after initialization: " << np << std::endl;

      for(int pindex = 0; pindex < np; ++pindex) {
          std::cout << "in particle loop after initialization" << std::endl;
          ParticleType& p = particles[pindex];
          const IntVect& iv = this->Index(p, lev);
          Real x, y, u, v, mFe, mFeO, mFe3O4, Hp, Tp, LFe, LFeO, LFe3O4;
          int phaseFe, phaseFeO, phaseFe3O4;
          x = p.pos(0);
          // y = p.pos(1);
          u = p.rdata(RealData::u);
          // v = p.rdata(RealData::v);
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

          std::cout << "got particle temperature post initialization" << std::endl;
      }
  }

  Redistribute();
}

void 
AmrLevelAdv::getParticleInfo(MultiFab& S_new)
{
  const int lev = 0;
  // const Geometry& geom = Geom(lev);
  // const Real* dx = geom.CellSize();
  // const Real* plo = geom.ProbLo();
  Real x, y, u, v, mFe, mFeO, mFe3O4, Hp, Tp, LFe, LFeO, LFe3O4, phaseFe, phaseFeO, phaseFe3O4;

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
        // y = p.pos(1);
        u = p.rdata(RealData::u);
        // v = p.rdata(RealData::v);
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

        std::cout << "x position is: " << x << ", x-velocity is: " << u << std::endl;

        std::cout << "got particle temperature in advance function" << std::endl;
    }
  }
  Redistribute();
}

void 
AmrLevelAdv::updateParticleInfo(MultiFab& Sborder, const double& dt, const double& dx, const double& dy)
{
  const int lev = 0;
  // const Geometry& geom = Geom(lev);
  // const Real* dx = geom.CellSize();
  // const Real* plo = geom.ProbLo();
  Real x, y, u, v, mFe, mFeO, mFe3O4, Hp, Tp, LFe, LFeO, LFe3O4;
  Real lx, ly;
  int  phaseFe, phaseFeO, phaseFe3O4;
  int  i=0, j=0, k=0;
  Vector<double> q(NUM_STATE);
  Vector<double> qSource(NUM_STATE,0);
  Vector<double> pReal(RealData::ncomps);
  Vector<double> pInt(IntData::ncomps);
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
        lx = (y - lo.x) / dx;
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
            pReal[i] = p.idata(h);
        }

        getSource(qSource,pSource,arr,i,j,k,pReal,pInt,dt,dx,dy,)


        // y = p.pos(1);
        // u = p.rdata(RealData::u);
        p.rdata(RealData::u) = p.rdata(RealData::u) + dt*10;
        u = p.rdata(RealData::u);
        // v = p.rdata(RealData::v);
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

        std::cout << "x position is: " << x << ", x-velocity is: " << u << std::endl;

        std::cout << "got particle temperature in advance function" << std::endl;
    }
  }
  Redistribute();
}
