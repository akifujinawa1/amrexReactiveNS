#ifndef EXACTFUNC_H_INCLUDED
#define EXACTFUNC_H_INCLUDED

// include dependencies here
#include <fstream> 
#include "AMReX_Vector.H"

using namespace amrex;

void updateExact(const int nSteps, const double dx, const int dim, std::ofstream &exact, Vector<Vector<double> >& arrExact);
double f(const double& p0,const double& pL,const double& pR, const double& aL, const double& aR, const double& vL,const double& vR,\
         const double& AL,const double& BL,const double& AR,const double& BR, const double& gamma);                 
double fL(const double& p0,const double& pL, const double& aL, const double& AL,const double& BL,const double& gamma);
double fR(const double& p0,const double& pR, const double& aR,const double& AR,const double& BR,const double& gamma);
double fp(const double& p0,const double& pL,const double& pR, const double& aL, const double& aR, const double& rhoL, const double& rhoR,\
         const double& AL,const double& BL,const double& AR,const double& BR, const double& gamma);

#endif