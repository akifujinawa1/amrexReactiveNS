#include "eulerFunc.H"
#include "AMReX_Vector.H"
#include "AMReX_REAL.H"

#include <math.h>
#include <iostream>
#include <algorithm> 

extern int enIC;
extern double Gamma;
extern int NUM_STATE;

namespace amrex;

// Use this file to write functions required for diffusive flux calculations

// grad(u) for viscous diffusion, stress tensor 
// grad(T) for heat diffusion
// grad(Y) for species mass diffusion
// tau for stress tensor


