#include "eulerFunc.H"
#include "diffusionFunc.H"
#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

#include <math.h>
#include <iostream>
#include <algorithm> 

extern int enIC;
extern double Gamma;
extern int NUM_STATE;
extern const int spacedim;
extern double R;
extern double M;
extern double Pr;

using namespace amrex;

// Use this file to write functions required for diffusive flux calculations

void updateViscous(MultiFab& Sborder, Array<MultiFab, SpaceDim>& fluxes, Vector<double> &qL, Vector<double> &qR, \
                   Vector<double> &qLlo, Vector<double> &qRlo, Vector<double> &qLhi, Vector<double> &qRhi, \
                   Vector<double> &viscSlice, const int &d, const double &dt, \
                   const double& dx, const double& dy, const int& SpaceDim, const int &viscous){
    switch (viscous)
    {
        case 0: 
        {}
        case 1: // Use central
        {
            const int iOffset = ( d == 0 ? 1 : 0);
            const int jOffset = ( d == 1 ? 1 : 0);
            const int kOffset = ( d == 2 ? 1 : 0);

            for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Dim3 lo = lbound(bx);
                const Dim3 hi = ubound(bx);

                // Indexable arrays for the data, and the directional flux
                // Based on the vertex-centred definition of the flux array, the
                // data array runs from e.g. [0,N] and the flux array from [0,N+1]
                const auto& arr = Sborder.array(mfi);
                const auto& fluxArrVisc = fluxes[d].array(mfi);
                const auto& fluxArrViscY = fluxes[d].array(mfi);

                // If slopelimiting has been enabled in the settings file, we run through the following set of loops
                // to compute boundary-reconstructed values -2023W2
                if (d == 0)
                {
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int j = lo.y; j <= hi.y; j++)
                        {
                            for(int i = lo.x; i <= hi.x+iOffset; i++) // at each i we calculate the left interface flux 
                            {           
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    qL[h]   = arr(i-1,j,k,h);   // cell to the left of cell i
                                    qR[h]   = arr(i  ,j,k,h);   // cell i 
                                    if (amrex::SpaceDim > 1)
                                    {
                                        qLlo[h] = arr(i-1,j-1,k,h);   // cell to the left lower diagonal of cell i
                                        qRlo[h] = arr(i  ,j-1,k,h);   // cell to the right lower diagonal of cell i 
                                        qLhi[h] = arr(i-1,j+1,k,h);   // cell to the left higher diagonal of cell i
                                        qRhi[h] = arr(i  ,j+1,k,h);   // cell to the right higher diagonal of cell i 
                                    }
                                }

                                if (amrex::SpaceDim == 1){ // x-direction viscous flux calculation in 1-D
                                    getViscFlux1D(viscSlice,qL,qR,dx);
                                    for(int h = 0; h < NUM_STATE; h++)
                                    {
                                    fluxArrVisc(i,j,k,h) = viscSlice[h]; // this is the viscous flux function for the left interface of cell i
                                    }
                                }
                                else { // x-direction viscous flux calculation if domain is 2-D
                                    getViscFlux2D(viscSlice,qL,qR,qLlo,qRlo,qLhi,qRhi,d,dx,dy);
                                    for(int h = 0; h < NUM_STATE; h++) {
                                    fluxArrVisc(i,j,k,h) = viscSlice[h];
                                    }
                                }
                            }
                        }
                    }
                }
                else // if solving for y-direction viscous fluxes (only applicable to 2-D)
                { 
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int i = lo.x; i <= hi.x; i++)
                        {
                            for(int j = lo.y; j <= hi.y+jOffset; j++)
                            {           
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    qL[h]   = arr(i,  j-1,k,h);   // cell below cell i,j
                                    qR[h]   = arr(i,  j  ,k,h);   // cell i,j 
                                    qLlo[h] = arr(i-1,j-1,k,h);   // cell to the left lower diagonal of cell i
                                    qRlo[h] = arr(i-1,j  ,k,h);   // cell to the right lower diagonal of cell i 
                                    qLhi[h] = arr(i+1,j-1,k,h);   // cell to the left higher diagonal of cell i
                                    qRhi[h] = arr(i+1,j  ,k,h);   // cell to the right higher diagonal of cell i 
                                }
                                getViscFlux2D(viscSlice,qL,qR,qLlo,qRlo,qLhi,qRhi,d,dx,dy);
                                for(int h = 0; h < NUM_STATE; h++) {
                                    fluxArrViscY(i,j,k,h) = viscSlice[h];
                                }
                            }
                        }
                    }
                }
                for(int k = lo.z; k <= hi.z; k++)
                {
                    for(int j = lo.y; j <= hi.y; j++)
                    {
                        for(int i = lo.x; i <= hi.x; i++)
                        {
                            // Conservative update formula
                            if (d == 0){ // x-direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) + (dt / dx) * (fluxArrVisc(i+iOffset, j, k,h) - fluxArrVisc(i,j,k,h));
                                }
                            }
                            else { // y=direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) + (dt / dy) * (fluxArrViscY(i, j+jOffset, k,h) - fluxArrViscY(i,j,k,h));
                                }
                            }
                           
                        }
                    }
                }
            } // close loop for patches
        } // close case 1 central
    } // close viscous switch case    

}

void getViscFlux1D(Vector<double>& viscSlice, const Vector<double>& qL,\
                   const Vector<double>& qR, const double& dx){

    double cV = 3.5*287;
    double rhoL, uL, enerL, epsL, TL, YL, pL, rhoR, uR, enerR, epsR, TR, YR, pR;
    double dudx, dTdx, dYdx, uAvg, TAvg, rhoAvg;
    double mu_avg, k_avg, D_avg;
    

    Vector<double> primL;
    Vector<double> primR;


    primL   = getPrim(qL);
    primR   = getPrim(qR);
    
    rhoL  = primL[0];
    uL    = primL[1];
    enerL = qL[3];
    pL    = primL[3];
    epsL  = specIntEner(rhoL, uL, 0, enerL);
    TL    = pL/((R/M)*rhoL);
    YL    = primL[4];
    
    rhoR  = primR[0];
    uR    = primR[1];
    enerR = qR[3];
    pR    = primR[3];
    epsR  = specIntEner(rhoR, uR, 0, enerR);
    TL    = pR/((R/M)*rhoR);
    YR    = primR[4];

    dudx = (uR-uL)/dx;
    dTdx = (TR-TL)/dx;
    dYdx = (YR-YL)/dx;
    uAvg = (uR+uL)/2.0;
    TAvg = (TR+TL)/2.0;
    rhoAvg = (rhoR+rhoL)/2.0;
    mu_avg = mu(TAvg,rhoAvg);
    k_avg  = k(TAvg,rhoAvg);
    D_avg  = D(TAvg,rhoAvg);

    viscSlice[0] = 0;
    viscSlice[1] = (4.0/3.0)*mu_avg*dudx;
    viscSlice[2] = 0;
    viscSlice[3] = uAvg*viscSlice[1] + k_avg*dTdx;
    viscSlice[4] = rhoAvg*D_avg*dYdx;

}

void getViscFlux2D(Vector<double>& viscSlice, const Vector<double>& qL, const Vector<double>& qR, \
                   const Vector<double>& qLlo, const Vector<double>& qRlo, const Vector<double>& qLhi, \
                   const Vector<double>& qRhi, const int& d, const double& dx, const double& dy){

    double cV = 3.5*287;
    double rhoL, uL, vL, enerL, epsL, TL, pL, YL, rhoR, uR, vR, enerR, epsR, TR, pR, YR;
    double rhoLlo, uLlo, vLlo, enerLlo, epsLlo, TLlo, pLlo, YLlo, rhoRlo, uRlo, vRlo, enerRlo, epsRlo, TRlo, pRlo, YRlo;
    double rhoLhi, uLhi, vLhi, enerLhi, epsLhi, TLhi, pLhi, YLhi, rhoRhi, uRhi, vRhi, enerRhi, epsRhi, TRhi, pRhi, YRhi;
    double dudx, dudy, dvdx, dvdy, dTdx, dTdy, dYdx, dYdy, uAvg, vAvg, TAvg, rhoAvg;
    // double mu = 3e-5;
    // double k  = 1.4;
    double mu_avg, k_avg, D_avg;

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
    enerL = qL[3];
    pL    = primL[3];
    epsL  = specIntEner(rhoL, uL, 0, enerL);
    TL    = pL/((R/M)*rhoL);
    YL    = primL[4];
    
    rhoR  = primR[0];
    uR    = primR[1];
    enerR = qR[3];
    pR    = primR[3];
    epsR  = specIntEner(rhoR, uR, 0, enerR);
    TL    = pR/((R/M)*rhoR);
    YR    = primR[4];

    rhoLlo  = primLlo[0];
    uLlo    = primLlo[1];
    enerLlo = qLlo[3];
    pLlo    = primLlo[3];
    epsLlo  = specIntEner(rhoLlo, uLlo, 0, enerLlo);
    TLlo    = pLlo/((R/M)*rhoLlo);
    YLlo    = primLlo[4];
    
    rhoRlo  = primRlo[0];
    uRlo    = primRlo[1];
    enerRlo = qRlo[3];
    pRlo    = primRlo[3];
    epsRlo  = specIntEner(rhoRlo, uRlo, 0, enerRlo);
    TLlo    = pRlo/((R/M)*rhoRlo);
    YRlo    = primRlo[4];

    rhoLhi  = primLhi[0];
    uLhi    = primLhi[1];
    enerLhi = qLhi[3];
    pLhi    = primLhi[3];
    epsLhi  = specIntEner(rhoLhi, uLhi, 0, enerLhi);
    TLhi    = pLhi/((R/M)*rhoLhi);
    YLhi    = primLhi[4];
    
    rhoRhi  = primRhi[0];
    uRhi    = primRhi[1];
    enerRhi = qRhi[3];
    pRhi    = primRhi[3];
    epsRhi  = specIntEner(rhoRhi, uRhi, 0, enerRhi);
    TLhi    = pRhi/((R/M)*rhoRhi);
    YRhi    = primRhi[4];


    if (d == 0) { // finding viscous fluxes for x-direction update

        dudx = (uR-uL)/dx;
        dvdx = (uR-uL)/dx;
        dTdx = (TR-TL)/dx;
        dYdx = (YR-YL)/dx;
        dudy = 0.5*((uRhi-uRlo)/(2*dy)+(uLhi-uLlo)/(2*dy));
        dvdy = 0.5*((vRhi-vRlo)/(2*dy)+(vLhi-vLlo)/(2*dy));
        uAvg = (uR+uL)/2.0;
        vAvg = (uR+uL)/2.0;
        TAvg = (TR+TL)/2.0;
        rhoAvg = (rhoR+rhoL)/2.0;
        mu_avg = mu(TAvg,rhoAvg);
        k_avg  = k(TAvg,rhoAvg);
        D_avg  = D(TAvg,rhoAvg);

        viscSlice[0] = 0;
        viscSlice[1] = 2.0*mu_avg*dudx - (2.0/3.0)*mu_avg*(dudx+dvdy);
        viscSlice[2] = mu_avg*(dvdx+dudy);
        viscSlice[3] = uAvg*viscSlice[1] + vAvg*viscSlice[2] + k_avg*dTdx;
        viscSlice[4] = rhoAvg*D_avg*dYdx;

    }

    else { // finding viscous fluxes for y-direction update

        dudy = (uR-uL)/dy;
        dvdy = (vR-vL)/dy;
        dTdy = (TR-TL)/dy;
        dYdy = (YR-YL)/dy;
        dudx = 0.5*((uRhi-uRlo)/(2*dx)+(uLhi-uLlo)/(2*dx));
        dvdx = 0.5*((vRhi-vRlo)/(2*dx)+(vLhi-vLlo)/(2*dx));
        uAvg = (uR+uL)/2.0;
        vAvg = (uR+uL)/2.0;
        TAvg = (TR+TL)/2.0;
        rhoAvg = (rhoR+rhoL)/2.0;
        mu_avg = mu(TAvg,rhoAvg);
        k_avg  = k(TAvg,rhoAvg);
        D_avg  = D(TAvg,rhoAvg);

        viscSlice[0] = 0;
        viscSlice[1] = mu_avg*(dvdx+dudy);
        viscSlice[2] = 2.0*mu_avg*dvdy - (2.0/3.0)*mu_avg*(dudx+dvdy);
        viscSlice[3] = uAvg*viscSlice[1] + vAvg*viscSlice[2] + k_avg*dTdy;
        viscSlice[4] = rhoAvg*D_avg*dYdy;


    }
    
}


double diffusiveSpeed(const Vector<double>& qL, const Vector<double>& qR){

    double rhoL, uL, enerL, epsL, TL, YL, pL, rhoR, uR, enerR, epsR, TR, YR, pR;
    double dudx, dTdx, dYdx, uAvg, TAvg, rhoAvg;
    double mu_avg, k_avg, D_avg;
    double maxLeft, maxRight;
    
    Vector<double> primL;
    Vector<double> primR;

    primL   = getPrim(qL);
    primR   = getPrim(qR);
    
    rhoL  = primL[0];
    uL    = primL[1];
    enerL = qL[3];
    pL    = primL[3];
    TL    = pL/((R/M)*rhoL);

    // std::cout << "temp: " << TL << std::endl;
    
    rhoR  = primR[0];
    uR    = primR[1];
    enerR = qR[3];
    pR    = primR[3];
    TL    = pR/((R/M)*rhoR);

    mu_avg = mu(TL,rhoL);
    k_avg  = k(TL,rhoL);
    D_avg  = D(TL,rhoL);

    maxLeft = 2*std::max(mu_avg/rhoL,mu_avg/(rhoL*Pr));

    mu_avg = mu(TR,rhoR);
    k_avg  = k(TR,rhoR);
    D_avg  = D(TR,rhoR);

    maxRight = 2*std::max(mu_avg/rhoR,mu_avg/(rhoR*Pr));

    // std::cout << "mu: " << mu_avg << ", max speed" << maxRight << std::endl;

    double maxSpeed = std::max(maxLeft,maxRight);
    return maxSpeed;

}


// functions to calculate transport properties via temperature fits go here

double mu(const double& T, const double& rho){
    double mu_val = (2.9*1e-4/rho)*pow(T,0.7);
    // std::cout << "mu: " << mu_val << std::endl;
    return mu_val;
}

double k(const double& T, const double& rho){
    double k_val = (2.9*1e-4/rho)*pow(T,0.7);
    return k_val;
}

double D(const double& T, const double& rho){
    double D_val = (2.9*1e-4/rho)*pow(T,0.7);
    return D_val;
}

