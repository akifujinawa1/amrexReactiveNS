#ifndef RECONSTRUCT_H_INCLUDED
#define RECONSTRUCT_H_INCLUDED

// include dependencies here
#include "AMReX_Vector.H"

extern int NUM_STATE;

using namespace amrex;

void reconstruct(Vector<double>& boundL, Vector<double>& boundR, const Vector<double>& slopeCells,\
                const Vector<Vector<double> >& u0, double w, double dx, double dt, const int d);

void getBoundsSLOPE(Vector<double>& boundL, Vector<double>& boundR, const Vector<double>& slopeCells,\
                    const Vector<Vector<double> >& u0, double w, const int d);

void localUpdate(Vector<double>& boundL, Vector<double>& boundR,\
                 double dx, double dt, const int d);

double getLimiter(double r, double rR);
double rval(double left, double center, double right);
double xsi_R(double r);

#endif