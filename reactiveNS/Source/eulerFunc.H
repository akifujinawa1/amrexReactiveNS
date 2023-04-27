#ifndef EULERFUNC_H
#define EULERFUNC_H

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"
#include "AMReX_Array.H"

#include <AMReX_MultiFab.H>

#include <vector>
#include <array>
#include <fstream> 


using namespace amrex;

extern const int spacedim;


void updateEuler(MultiFab& Sborder, Array<MultiFab, SpaceDim>& fluxes, Vector<double> &qL, Vector<double> &qR, \
                 Vector<double> &fluxvals, const int &SpaceDim, const double &dt, \
                 const double& dx, const double& dy, const int &euler);

Vector<double> getEulerFlux(const Vector<double>& q, const int& dir);
Vector<double> getPrim(const Vector<double>& q);
Vector<double> getCons(const Vector<double>& prim);
Vector<double> setIC(const int dim);
Vector<double> HLLCflux(const Vector<double>& qL, const Vector<double>& qR, const int& dir);
Vector<double> HLLCstarFlux(const Vector<double>& f, const Vector<double>& q,\
                                 const double& s, const double& sStar, const int& dir);

void wavespeedEstimate(const Vector<double>& qL, const Vector<double>& qR,\
                       double& sL, double& sR, double& sStar, const int dir);
void getStopTime(int enIC, amrex::Real& stop_time);


//double timeStepRP(const Vector<double>& qL, const Vector<double>& qR, const int& counter);
double soundSpeed(const double& p, const double& rho);
double energy(double rho, double vx, double vy, double p);
double pressure(double rho, double vx, double vy, double energy);
double specIntEner(double rho, double vx, double vy, double energy);

// enum IC {
//     toro1    = 0,
//     toro2    = 1,
//     toro3    = 2,
//     toro4    = 3,
//     toro5    = 4, 
//     cylExp   = 5,
//     tororot  = 6
// };

enum Limiter {
    Superbee    = 0,
    vanLeer     = 1,
    vanAlbada   = 2,
    minbee      = 3
};


#endif