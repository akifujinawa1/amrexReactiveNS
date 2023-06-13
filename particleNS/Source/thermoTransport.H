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
double cpMixFe(const double& cpO2, const double& cpN2, const double& cpFe, const double& cpFeO, \
               const double& Y_O2, const double& Y_N2, const double& Y_Fe, const double& Y_FeO);
double kO2(const double& Tg);
double kN2(const double& Tg);
double kFe(const double& Tg);
double kFeO(const double& Tg);
double kMix(const double& kO2, const double& kN2, const double& kFe, const double& kFeO,\
            const double& Y_O2, const double& Y_N2, const double& Y_Fe, const double& Y_FeO);
double muO2(const double& Tg);
double muN2(const double& Tg);
double muFe(const double& Tg);
double muFeO(const double& Tg);
double muMix(const double& muO2, const double& muN2, const double& muFe, const double& muFeO, \
             const double& YO2, const double& YN2, const double& YFe, const double& YFeO);
double Hgas(const double& rho, const double& YO2, const double& YN2, const double& Tg);
double Tg(const double& rho, const double& u, const double& v, \
          const double& YO2, const double& YN2, const double& ener);
double getMavg(const double& YO2, const double& YN2);
double omega(const double& T, const double& eps);
Vector<double> getMoleFractions(const double& Y_O2, const double& Y_N2, const double& Y_Fe, const double& Y_FeO);
Vector<double> getMassFractions(const double& X_O2, const double& X_N2, const double& X_Fe, const double& X_FeO);
Vector<double> getMixDiffCoeffs(const double& T, const double& p, const double& YO2, \
                                const double& YN2, const double& YFe, const double& YFeO);


#endif