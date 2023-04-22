#ifndef DIFFUSIONFUNC_H
#define DIFFUSIONFUNC_H

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"

#include <vector>
#include <array>
#include <fstream> 

using namespace amrex;

void getViscFlux1D(Vector<double>& viscSlice, const Vector<double>& qL,\
                   const Vector<double>& qR, const double& dx);
void getViscFlux2D(Vector<double>& viscSlice, const Vector<double>& qL, const Vector<double>& qR, \
                   const Vector<double>& qLlo, const Vector<double>& qRlo, const Vector<double>& qLhi, \
                   const Vector<double>& qRhi, const int& d, const double& dx, const double& dy);
double mu(const double& T);
double k(const double& T);

#endif