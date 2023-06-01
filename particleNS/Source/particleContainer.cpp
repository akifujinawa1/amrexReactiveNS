#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

// #include "particleContainer.H"
#include "AmrLevelAdv.H"
#include "particleFunc.H"


using namespace amrex;

// namespace {

//     void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
//     {
//         int nx = nppc[0];
//         int ny = nppc[1];
//         int nz = nppc[2];

//         int ix_part = i_part/(ny * nz);
//         int iy_part = (i_part % (ny * nz)) % ny;
//         int iz_part = (i_part % (ny * nz)) / ny;

//         r[0] = (0.5+ix_part)/nx;
//         r[1] = (0.5+iy_part)/ny;
//         r[2] = (0.5+iz_part)/nz;
//     }

//     void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std) {
//         Real ux_th = amrex::RandomNormal(0.0, u_std);
//         Real uy_th = amrex::RandomNormal(0.0, u_std);
//         Real uz_th = amrex::RandomNormal(0.0, u_std);

//         u[0] = u_mean + ux_th;
//         u[1] = u_mean + uy_th;
//         u[2] = u_mean + uz_th;
//     }
// }

// Empty constructor, builds invalid object to pass to initData
// SingleParticleContainer::SingleParticleContainer ()
// {
// }

// Basic constructor

SingleParticleContainer::SingleParticleContainer(const Geometry            & a_geom,
                                     const DistributionMapping & a_dmap,
                                     const BoxArray            & a_ba)
    : ParticleContainer<RealData::ncomps, IntData::ncomps> (a_geom, a_dmap, a_ba)
{
}

void SingleParticleContainer::initParticles () {

    // for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    const int lev = 0;
    using MyParIter = amrex::ParIter<RealData::ncomps, IntData::ncomps>;
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        // const Box& tile_box  = mfi.tilebox();
        // const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        // const int grid_id = mfi.index();
        // const int tile_id = mfi.LocalTileIndex();
        // auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)]

        auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];

        ParticleType p;
        p.id()   = ParticleType::NextID();
        p.cpu()  = ParallelDescriptor::MyProc();

        p.pos(0) = 0.5;
        // p.pos(1) =  0.0;

        double mFe0,mFeO0,mFe3O40,energy0;
        particleInit(mFe0,mFeO0,mFe3O40,energy0);

        p.rdata(RealData::u)     = 0;         // Store the x-velocity here (unused if simulation is 1-D)
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



// void SingleParticleContainer::getParticleInfo () {
        
//     std::cout << "in getParticleInfo" << std::endl;

//     const int lev = 0;
//     std::cout << "accessing Geom(lev)" << std::endl;
//     const Geometry& geom = Geom(lev);
//     // const Real* dx = geom.CellSize();
//     // const Real* plo = geom.ProbLo();
//     Real x, u, mFe, mFeO, mFe3O4, Hp, Tp, LFe, LFeO, LFe3O4, phaseFe, phaseFeO, phaseFe3O4;

//     const int grid_id = 0;
//     const int tile_id = 0;

//     auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
//     auto& particles = particle_tile.GetArrayOfStructs();
//     ParticleType& p = particles[0];
//     std::cout << "accessing 0th particle" << std::endl;

//     x = p.pos(0);
//     std::cout << "accessed 0th particle position" << std::endl;
//     // y = p.pos(1);
//     u = p.rdata(RealData::u);
//     // v = p.rdata(RealData::v);
//     mFe    = p.rdata(RealData::mFe);
//     mFeO   = p.rdata(RealData::mFeO);
//     mFe3O4 = p.rdata(RealData::mFe3O4);
//     Hp     = p.rdata(RealData::Hp);
//     LFe    = p.rdata(RealData::LFe);
//     LFeO   = p.rdata(RealData::LFeO);
//     LFe3O4 = p.rdata(RealData::LFe3O4);

//     phaseFe    = p.idata(IntData::Fe);
//     phaseFeO   = p.idata(IntData::FeO);
//     phaseFe3O4 = p.idata(IntData::Fe3O4);

//     Tp = Tparticle(mFe,mFeO,mFe3O4,Hp,phaseFe,phaseFeO,phaseFe3O4,LFe,LFeO,LFe3O4);

//     std::cout << "got particle temperature" << std::endl;
// }


// void SingleParticleContainer::InitParticles(MultiFab& S_new, const int num_particles){

//     const int lev = 0;
//     const Geometry& geom = Geom(lev);
//     const Real* dx = geom.CellSize();
//     const Real* plo = geom.ProbLo();

//     const int num_ppc = num_particles;

//     // const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
//     //                                  *a_num_particles_per_cell[1],
//     //                                  *a_num_particles_per_cell[2]);

//     for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
//     // for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
//     {
//         const Box& tile_box  = mfi.tilebox();
//         const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
//         const int grid_id = mfi.index();
//         const int tile_id = mfi.LocalTileIndex();
//         auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

//         // for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
//         // {
//             // for (int i_part=0; i_part<num_ppc;i_part++) {
//                 // Real r[3];
//                 // Real u[3];

//                 // get_position_unit_cell(r, a_num_particles_per_cell, i_part);
//                 // get_gaussian_random_momentum(u, 0.0, 0.1);

//                 // Real x = plo[0] + (iv[0] + r[0])*dx[0];
//                 // Real y = plo[1] + (iv[1] + r[1])*dx[1];
//                 // Real z = plo[2] + (iv[2] + r[2])*dx[2];

//                 Real x = plo[0] + 0.5*dx[0];
//                 Real y = 0;
//                 if (amrex::SpaceDim == 2){
//                     Real y = plo[1] + 0.5*dx[1];
//                 }
                
//                 ParticleType p;
//                 p.id()  = ParticleType::NextID();
//                 p.cpu() = ParallelDescriptor::MyProc();

//                 std::cout << "particle ID: " << p.id() << std::endl;

//                 p.pos(0) = x;
//                 p.pos(1) = y;

//                 double mFe0,mFeO0,mFe3O40,energy0;
//                 particleInit(mFe0,mFeO0,mFe3O40,energy0);

//                 p.rdata(RealData::u)     = 0;         // Store the x-velocity here (unused if simulation is 1-D)
//                 p.rdata(RealData::v)     = 0.0;       // Store the y-velocity here (unused if simulation is 1-D)
//                 p.rdata(RealData::w)     = 0.0;       // Store the z-velocity here (unused if simulation is 1-D)
//                 p.rdata(RealData::mFe)    = mFe0;      // Store the initial Fe mass here
//                 p.rdata(RealData::mFeO)   = mFeO0;     // Store the initial FeO mass here
//                 p.rdata(RealData::mFe3O4) = mFe3O40;   // Store the initial Fe3O4 mass here
//                 p.rdata(RealData::Hp)     = energy0;   // Store the initial particle enthalpy here
//                 p.rdata(RealData::LFe)    = 0.0;       // Store the Fe melt progress variable here
//                 p.rdata(RealData::LFeO)   = 0.0;       // Store the FeO melt progress variable here
//                 p.rdata(RealData::LFe3O4) = 0.0;       // Store the Fe3O4 melt progress variable here

//                 p.idata(IntData::Fe)     = 0;         // Store the Fe melt flag variable here
//                 p.idata(IntData::FeO)    = 0;         // Store the FeO melt flag variable here
//                 p.idata(IntData::Fe3O4)  = 0;         // Store the Fe3O4 melt flag variable here

//                 AMREX_ASSERT(this->Index(p, lev) == iv);

//                 // std::pair<int,int> key {grid_id,tile_id};
//                 // auto& particle_tile = GetParticles(0)[key];

//                 particle_tile.push_back(p);
//                 // particle_tile.push_back_real(attribs);

//                 // particle_tile.push_back(p);
//                 std::cout << "grid id: " << grid_id << ", tile id: " << tile_id << std::endl;
//             // }
//         // }
//         // auto& particles = GetParticles(lev);
//         // auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
//         // auto old_size = particle_tile.GetArrayOfStructs().size();
//         // auto new_size = old_size + host_particles.size();
//         // particle_tile.resize(new_size);
//     }
//     // Redistribute();
// }



void SingleParticleContainer::getParticleInfo(MultiFab& S_new){

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    Real x, y, u, v, mFe, mFeO, mFe3O4, Hp, Tp, LFe, LFeO, LFe3O4, phaseFe, phaseFeO, phaseFe3O4;

    std::cout << "in getParticleInfo" << std::endl;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();
        std::cout << "number of particles in grid after initialization: " << np << std::endl;


    // for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    // {
    //     std::cout << "in MyParIter for loop" << std::endl;
    //     const int grid_id = pti.index();
    //     const int tile_id = pti.LocalTileIndex();
    // //     // const int grid_id = 0;
    // //     // const int tile_id = 0;
    //     auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    //     auto& particles = particle_tile.GetArrayOfStructs();
    //     const int np = particles.numParticles();

    //     std::cout << "grid_id, tile_id, np: " << grid_id << " " << tile_id << " " << np << std::endl;
    //     // np = 1;
    //     // for (const auto& p : particles) {
    //     for(int pindex = 0; pindex < np; ++pindex) {
    //         ParticleType& p = particles[pindex];
    //         const IntVect& iv = this->Index(p, lev);

    // for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){

    //     // Index of grid (= box)
    //     const int grid_id = mfi.index();
    //     // Index of tile within the grid
    //     const int tile_id = mfi.LocalTileIndex();

    //     auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

    //     auto& particles = particle_tile.GetArrayOfStructs();
    //     const int np = particles.numParticles();


        // std::cout << "number of particles: " << np << std::endl;
        // Get GPU-friendly arrays of particle data
        // auto& ptile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        // ParticleType* p = particle_tile.GetArrayOfStructs()().data();
        // Only need attribs (i.e., SoA data)
        // auto& soa = ptile.GetStructOfArrays();
        // As an example, let's get the ux momentum
        // const ParticleReal * const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();

    // for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
    //     std::cout << "in MyParIter" << std::endl;
    //     auto& particles = pti.GetArrayOfStructs();
        // const int np    = pti.numParticles();
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

            std::cout << "got particle temperature" << std::endl;
        }
    }
}







// void
// CellSortedParticleContainer::UpdateCellVectors()
// {
//     BL_PROFILE("CellSortedParticleContainer::UpdateCellVectors");

//     const int lev = 0;

//     bool needs_update = false;
//     if (! m_vectors_initialized)
//     {
//         // this is the first call, so we must update
//         m_vectors_initialized = true;
//         needs_update = true;
//     }
//     else if ((m_BARef != this->ParticleBoxArray(lev).getRefID()) ||
//              (m_DMRef != this->ParticleDistributionMap(lev).getRefID()))
//     {
//         // the grids have changed, so we must update
//         m_BARef = this->ParticleBoxArray(lev).getRefID();
//         m_DMRef = this->ParticleDistributionMap(lev).getRefID();
//         needs_update = true;
//     }

//     if (! needs_update) return;

//     // clear old data
//     m_cell_vectors.clear();
//     m_vector_size.clear();
//     m_vector_ptrs.clear();

//     // allocate storage for cell vectors. NOTE - do not tile this loop
//     for(MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi)
//     {
//         const Box& box = mfi.validbox();
//         const int grid_id = mfi.index();
//         m_cell_vectors[grid_id].resize(box);
//         m_vector_size[grid_id].resize(box);
//         m_vector_ptrs[grid_id].resize(box);
//     }

//     // insert particles into vectors - this can be tiled
// #ifdef AMREX_USE_OMP
// #pragma omp parallel
// #endif
//     for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
//     {
//         auto& particles = pti.GetArrayOfStructs();
//         const int np    = pti.numParticles();
//         for(int pindex = 0; pindex < np; ++pindex) {
//             ParticleType& p = particles[pindex];
//             const IntVect& iv = this->Index(p, lev);
//             p.idata(IntData::sorted) = 1;
//             p.idata(IntData::i) = iv[0];
//             p.idata(IntData::j) = iv[1];
//             p.idata(IntData::k) = iv[2];
//             // note - use 1-based indexing for convenience with Fortran
//             m_cell_vectors[pti.index()](iv).push_back(pindex + 1);
//         }
//     }

//     UpdateFortranStructures();
// }

// void
// CellSortedParticleContainer::UpdateFortranStructures()
// {
//     BL_PROFILE("CellSortedParticleContainer::UpdateFortranStructures");

//     const int lev = 0;

// #ifdef AMREX_USE_OMP
// #pragma omp parallel
// #endif
//     for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
//     {
//         const Box& tile_box  = mfi.tilebox();
//         const int grid_id = mfi.index();
//         for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
//         {
//             m_vector_size[grid_id](iv) = m_cell_vectors[grid_id](iv).size();
//             m_vector_ptrs[grid_id](iv) = m_cell_vectors[grid_id](iv).data();
//         }
//     }
// }

// void
// CellSortedParticleContainer::MoveParticles()
// {
//     BL_PROFILE("CellSortedParticleContainer::MoveParticles()");

//     UpdateCellVectors();

//     const int lev = 0;
//     const Real* dx = Geom(lev).CellSize();
//     const Real* plo = Geom(lev).ProbLo();
//     const Real dt = 0.1;

// #ifdef AMREX_USE_OMP
// #pragma omp parallel
// #endif
//     for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
//     {
//         const int grid_id = pti.index();
//         const int tile_id = pti.LocalTileIndex();
//         const Box& tile_box  = pti.tilebox();

//         auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
//         auto& particles = particle_tile.GetArrayOfStructs();
//         const int np = particles.numParticles();

//         move_particles(particles.data(), &np,
//                        tile_box.loVect(), tile_box.hiVect(),
//                        m_vector_ptrs[grid_id].dataPtr(),
//                        m_vector_size[grid_id].dataPtr(),
//                        m_vector_ptrs[grid_id].loVect(),
//                        m_vector_ptrs[grid_id].hiVect(),
//                        plo, dx, &dt);

//         // resize particle vectors after call to move_particles
//         for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
//         {
//             const auto new_size = m_vector_size[grid_id](iv);
//             auto& pvec = m_cell_vectors[grid_id](iv);
//             pvec.resize(new_size);
//         }
//     }
// }

// void
// CellSortedParticleContainer::ReBin()
// {
//     BL_PROFILE("CellSortedParticleContainer::ReBin()");

//     const int lev = 0;

// #ifdef AMREX_USE_OMP
// #pragma omp parallel
// #endif
//     for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
//     {
//         const int grid_id = pti.index();
//         const int tile_id = pti.LocalTileIndex();

//         auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
//         auto& particles = particle_tile.GetArrayOfStructs();
//         const int np = particles.numParticles();
//         for(int pindex = 0; pindex < np; ++pindex)
//         {
//             ParticleType& p = particles[pindex];
//             if (p.idata(IntData::sorted)) continue;
//             const IntVect& iv = this->Index(p, lev);
//             p.idata(IntData::sorted) = 1;
//             p.idata(IntData::i) = iv[0];
//             p.idata(IntData::j) = iv[1];
//             p.idata(IntData::k) = iv[2];
//             // note - use 1-based indexing for convenience with Fortran
//             m_cell_vectors[pti.index()](iv).push_back(pindex + 1);
//         }
//     }

//     UpdateFortranStructures();
// }
