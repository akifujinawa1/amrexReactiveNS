#ifndef EULERFUNC_H
#define EULERFUNC_H

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"

#include <vector>
#include <array>
#include <fstream> 


amrex::Vector<double> getEulerFlux(const amrex::Vector<double>& q, const int& dir);
amrex::Vector<double> getPrim(const amrex::Vector<double>& q);
amrex::Vector<double> getCons(const amrex::Vector<double>& prim);
amrex::Vector<double> setIC(const int dim);
amrex::Vector<double> HLLCflux(const amrex::Vector<double>& qL, const amrex::Vector<double>& qR, const int& dir);
amrex::Vector<double> HLLCstarFlux(const amrex::Vector<double>& f, const amrex::Vector<double>& q,\
                                 const double& s, const double& sStar, const int& dir);

void wavespeedEstimate(const amrex::Vector<double>& qL, const amrex::Vector<double>& qR,\
                       double& sL, double& sR, double& sStar, const int dir);
void getStopTime(int enIC, amrex::Real& stop_time);

//double timeStepRP(const amrex::Vector<double>& qL, const amrex::Vector<double>& qR, const int& counter);
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