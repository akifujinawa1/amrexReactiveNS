#include "eulerFunc.H"
#include "diffusionFunc.H"
#include "AMReX_Vector.H"
#include "AMReX_REAL.H"

#include <math.h>
#include <iostream>
#include <algorithm> 

extern int enIC;
extern double Gamma;
extern int NUM_STATE;

using namespace amrex;

// Use this file to write functions required for diffusive flux calculations

void getViscFlux1D(Vector<double>& viscSlice, const Vector<double>& qL,\
                   const Vector<double>& qR, const double& dx){

    double cV = 3.5*287;
    double rhoL, uL, enerL, epsL, TL, rhoR, uR, enerR, epsR, TR;
    double dudx, dTdx, uAvg, TAvg;
    double mu_avg, k_avg;
    // double mu = 3e-5;
    // double k  = 1.4;
    // double mu, k;

    Vector<double> primL;
    Vector<double> primR;


    primL   = getPrim(qL);
    primR   = getPrim(qR);
    
    rhoL  = primL[0];
    uL    = primL[1];
    enerL = primL[3];
    epsL  = specIntEner(rhoL, uL, 0, enerL);
    TL    = epsL/cV;
    
    rhoR  = primR[0];
    uR    = primR[1];
    enerR = primR[3];
    epsR  = specIntEner(rhoR, uR, 0, enerR);
    TR    = epsR/cV;

    dudx = (uR-uL)/dx;
    dTdx = (TR-TL)/dx;
    uAvg = (uR+uL)/2.0;
    TAvg = (TR+TL)/2.0;
    mu_avg = mu(TAvg);
    k_avg  = k(TAvg);

    viscSlice[0] = 0;
    viscSlice[1] = (4.0/3.0)*mu_avg*dudx;
    viscSlice[2] = 0;
    viscSlice[3] = uAvg*viscSlice[1] + k_avg*dTdx;

}

void getViscFlux2D(Vector<double>& viscSlice, const Vector<double>& qL, const Vector<double>& qR, \
                   const Vector<double>& qLlo, const Vector<double>& qRlo, const Vector<double>& qLhi, \
                   const Vector<double>& qRhi, const int& d, const double& dx, const double& dy){

    double cV = 3.5*287;
    double rhoL, uL, vL, enerL, epsL, TL, rhoR, uR, vR, enerR, epsR, TR;
    double rhoLlo, uLlo, vLlo, enerLlo, epsLlo, TLlo, rhoRlo, uRlo, vRlo, enerRlo, epsRlo, TRlo;
    double rhoLhi, uLhi, vLhi, enerLhi, epsLhi, TLhi, rhoRhi, uRhi, vRhi, enerRhi, epsRhi, TRhi;
    double dudx, dudy, dvdx, dvdy, dTdx, dTdy, uAvg, vAvg, TAvg;
    // double mu = 3e-5;
    // double k  = 1.4;
    double mu_avg, k_avg;

    Vector<double> primL;
    Vector<double> primR;
    Vector<double> primLlo;
    Vector<double> primRlo;
    Vector<double> primLhi;
    Vector<double> primRhi;

    primL   = getPrim(qL);
    primR   = getPrim(qR);
    primLlo = getPrim(qLlo);
    primRlo = getPrim(qRlo);
    primLhi = getPrim(qLhi);
    primRhi = getPrim(qRhi);
    
    rhoL  = primL[0];
    uL    = primL[1];
    enerL = primL[3];
    epsL  = specIntEner(rhoL, uL, 0, enerL);
    TL    = epsL/cV;
    
    rhoR  = primR[0];
    uR    = primR[1];
    enerR = primR[3];
    epsR  = specIntEner(rhoR, uR, 0, enerR);
    TR    = epsR/cV;

    rhoLlo  = primLlo[0];
    uLlo    = primLlo[1];
    enerLlo = primLlo[3];
    epsLlo  = specIntEner(rhoLlo, uLlo, 0, enerLlo);
    TLlo    = epsLlo/cV;
    
    rhoRlo  = primRlo[0];
    uRlo    = primRlo[1];
    enerRlo = primRlo[3];
    epsRlo  = specIntEner(rhoRlo, uRlo, 0, enerRlo);
    TRlo    = epsRlo/cV;

    rhoLhi  = primLhi[0];
    uLhi    = primLhi[1];
    enerLhi = primLhi[3];
    epsLhi  = specIntEner(rhoLhi, uLhi, 0, enerLhi);
    TLhi    = epsLhi/cV;
    
    rhoRhi  = primRhi[0];
    uRhi    = primRhi[1];
    enerRhi = primRhi[3];
    epsRhi  = specIntEner(rhoRhi, uRhi, 0, enerRhi);
    TRhi    = epsRhi/cV;

    if (d == 0) { // finding viscous fluxes for x-direction update

        dudx = (uR-uL)/dx;
        dvdx = (uR-uL)/dx;
        dTdx = (TR-TL)/dx;
        dudy = 0.5*((uRhi-uRlo)/(2*dy)+(uLhi-uLlo)/(2*dy));
        dvdy = 0.5*((vRhi-vRlo)/(2*dy)+(vLhi-vLlo)/(2*dy));
        uAvg = (uR+uL)/2.0;
        vAvg = (uR+uL)/2.0;
        TAvg = (TR+TL)/2.0;
        mu_avg = mu(TAvg);
        k_avg  = k(TAvg);

        viscSlice[0] = 0;
        viscSlice[1] = 2.0*mu_avg*dudx - (2.0/3.0)*mu_avg*(dudx+dvdy);
        viscSlice[2] = mu_avg*(dvdx+dudy);
        viscSlice[3] = uAvg*viscSlice[1] + vAvg*viscSlice[2] + k_avg*dTdx;

    }

    else { // finding viscous fluxes for y-direction update

        dudy = (uR-uL)/dy;
        dvdy = (vR-vL)/dy;
        dTdy = (TR-TL)/dy;
        dudx = 0.5*((uRhi-uRlo)/(2*dx)+(uLhi-uLlo)/(2*dx));
        dvdx = 0.5*((vRhi-vRlo)/(2*dx)+(vLhi-vLlo)/(2*dx));
        uAvg = (uR+uL)/2.0;
        vAvg = (uR+uL)/2.0;
        TAvg = (TR+TL)/2.0;
        mu_avg = mu(TAvg);
        k_avg  = k(TAvg);

        viscSlice[0] = 0;
        viscSlice[1] = mu_avg*(dvdx+dudy);
        viscSlice[2] = 2.0*mu_avg*dvdy - (2.0/3.0)*mu_avg*(dudx+dvdy);
        viscSlice[3] = uAvg*viscSlice[1] + vAvg*viscSlice[2] + k_avg*dTdy;

    }
    
}


// functions to calculate transport properties via temperature fits go here

double mu(const double& T){
    double mu_val = 3e-5;
    return mu_val;
}

double k(const double& T){
    double k_val = 1.4;
    return k_val;
}

// double mu = 3e-5;
    // double k  = 1.4;

// primLlo = getPrim(qLlo);
    // primRlo = getPrim(qRlo);
    // primLhi = getPrim(qLhi);
    // primRhi = getPrim(qRhi);