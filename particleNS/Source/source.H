// This file contains various function headers that are used to compute the source terms
// for the augmented Euler equations. Used with source.cpp

#ifndef SOURCE_H_INCLUDED
#define SOURCE_H_INCLUDED

// include dependencies here
#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

extern int NUM_STATE;

using namespace amrex;

void updateSource(MultiFab& Sborder, Vector<double>& q, Vector<double>& sourceVec, const double &dt, const int &source);
void getSource(const Vector<double>& q, Vector<double>& sourceVec, const int& source, const double& dt);
void EulerEx(Vector<double>& source, const Vector<double>& q, const double& dt);
void rk4(Vector<double>& source, const Vector<double>& q, const double& dt);
double getEnergySource(const double& rho, const double& lambda, const double& T);
double getReactiveSource(const double& rho, const double& lambda, const double& T);


#endif