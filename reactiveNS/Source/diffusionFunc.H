#ifndef DIFFUSIONFUNC_H
#define DIFFUSIONFUNC_H

#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"

#include <vector>
#include <array>
#include <fstream> 

extern const int spacedim;

using namespace amrex;

void updateViscous(MultiFab& Sborder, Array<MultiFab, SpaceDim>& fluxes, Vector<double> &qL, Vector<double> &qR, \
                   Vector<double> &qLlo, Vector<double> &qRlo, Vector<double> &qLhi, Vector<double> &qRhi, \
                   Vector<double> &viscSlice, Vector<double> &fluxvals, const int &d, const double &dt, \
                   const double& dx, const double& dy, const int& SpaceDim, const int &viscous);
                   
void getViscFlux1D(Vector<double>& viscSlice, const Vector<double>& qL,\
                   const Vector<double>& qR, const double& dx);
void getViscFlux2D(Vector<double>& viscSlice, const Vector<double>& qL, const Vector<double>& qR, \
                   const Vector<double>& qLlo, const Vector<double>& qRlo, const Vector<double>& qLhi, \
                   const Vector<double>& qRhi, const int& d, const double& dx, const double& dy);
double mu(const double& T);
double k(const double& T);

#endif