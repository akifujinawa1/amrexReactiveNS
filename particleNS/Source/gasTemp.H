#ifndef GASTEMP_H
#define GASTEMP_H


double hO2(const double& Tg);
double hN2(const double& Tg);
double cpO2(const double& Tg);
double cpN2(const double& Tg);
double cpMix(const double& cpO2, const double& cpN2, const double& Y_O2, const double& Y_N2);
double Hgas(const double& rho, const double& YO2, const double& YN2, const double& Tg);
double Tg(const double& rho, const double& u, const double& v, \
          const double& YO2, const double& YN2, const double& ener);
double getMavg(const double& YO2, const double& YN2);


#endif