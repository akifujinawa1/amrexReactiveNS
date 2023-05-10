#ifndef THERMOTRANSPORT_H
#define THERMOTRANSPORT_H

#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"

#include <vector>

using namespace amrex;

double hO2(const double& Tg);
double hN2(const double& Tg);
double cpO2(const double& Tg);
double cpN2(const double& Tg);
double cpMix(const double& cpO2, const double& cpN2, const double& Y_O2, const double& Y_N2);
double kO2(const double& Tg);
double kN2(const double& Tg);
double kMix(const double& kO2, const double& kN2, const double& Y_O2, const double& Y_N2);
double muO2(const double& Tg);
double muN2(const double& Tg);
double muMix(const double& muO2, const double& muN2, const double& YO2, const double& YN2);
double Hgas(const double& rho, const double& YO2, const double& YN2, const double& Tg);
double Tg(const double& rho, const double& u, const double& v, \
          const double& YO2, const double& YN2, const double& ener);
double getMavg(const double& YO2, const double& YN2);
Vector<double> getMoleFractions(const double& YO2, const double& YN2);
Vector<double> getMassFractions(const double& XO2, const double& XN2);


#endif