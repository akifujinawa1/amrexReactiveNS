#ifndef PARTICLE_CONTAINER_H_
#define PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

struct RealData
{
    enum {
        u = 0,
        v,
        w,
        mFe,
        mFeO,
        mFe3O4,
        Hp,
        LFe,
        LFeO,
        LFe3O4,
        ncomps
    };
};

struct IntData
{
    enum {
        Fe,
        FeO,
        Fe3O4,
        ncomps
    };
};

// class MyParIter
//     : public amrex::ParIter<RealData::ncomps, IntData::ncomps>
// {
// public:
//     using amrex::ParIter<RealData::ncomps, IntData::ncomps>::ParIter;
// };

class SingleParticleContainer
    : public amrex::ParticleContainer<RealData::ncomps, IntData::ncomps>
{
public:

    SingleParticleContainer (const amrex::Geometry            & a_geom,
                             const amrex::DistributionMapping & a_dmap,
                             const amrex::BoxArray            & a_ba);

    using MyParIter = amrex::ParIter<RealData::ncomps, IntData::ncomps>;

    static void initParticles();  //amrex::MultiFab& S_new, const int num_particles

    void getParticleInfo(amrex::MultiFab& S_new);  //

protected:

private:

};

#endif