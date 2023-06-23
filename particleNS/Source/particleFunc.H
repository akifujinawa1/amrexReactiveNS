#ifndef PARTICLEFUNC_H
#define PARTICLEFUNC_H

#include "AMReX_Vector.H"

void particleInit(double& mFe0, double& mFeO0, double& mFe3O40, double& energy0);
double Cdrag(double Re);
double Hparticle(const double& mFe, const double& mFeO, const double& mFe3O4, const double& Tp, \
                 const int& phaseFe, const int& phaseFeO, const int& phaseFe3O4);
double hFeS(const double& Tp);
double hFeL(const double& Tp);
double hFeG(const double& Tp);
double hFeOS(const double& Tp);
double hFeOL(const double& Tp);
double hFeOG(const double& Tp);
double hFe3O4S(const double& Tp);
double hFe3O4L(const double& Tp);
double cpparticle(const double& mFe, const double& mFeO, const double& mFe3O4, const double& Tp, \
                  const int& phaseFe, const int& phaseFeO, const int& phaseFe3O4);
double cpFeS(const double& Tp);
double cpFeL(const double& Tp);
double cpFeG(const double& Tp);
double cpFeOS(const double& Tp);
double cpFeOL(const double& Tp);
double cpFeOG(const double& Tp);
double cpFe3O4S(const double& Tp);
double cpFe3O4L(const double& Tp);
double Tparticle(const double& mFe, const double& mFeO, const double& mFe3O4, const double& Hp, \
                 int& phaseFe, int& phaseFeO, int& phaseFe3O4, \
                 const double& LFe, const double& LFeO, const double& LFe3O4);
void getSource(Vector<double>& qSource, Vector<double>& pSource, const Vector<double>& q, Vector<double>& pReal, \
               Vector<int>& pInt, const double& dt, const double& dx, const double& dy);

Vector<double> getOxidationRates(const double& mFe, const double& mFeO, const double& mFe3O4, \
                                 const double& Tp,  const double& rp);
double filmAverage(const double& Tp, const double& Tg);

#endif