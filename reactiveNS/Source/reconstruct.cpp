// This file contains various functions that are used to compute the linear boundary
// reconstructed values for the MUSCL-Hancock method.


// include header files here

#include "eulerFunc.H"
#include "recon.H"
#include "AMReX_Vector.H"

// include dependencies here

#include <math.h>
#include <iostream>
#include <algorithm> 

extern int enLimiter;
extern int NUM_STATE;

using namespace amrex;


void reconstruct(Vector<double>& boundL, Vector<double>& boundR, const Vector<double>& slopeCells,\
                const Vector<Vector<double> >& u0, double w, double dx, double dt, const int d){
    
    // This is the main function that is called from AmrLevelAdv.cpp. Updates the vectors
    // 'boundL' and 'boundR' by reference by getting the slope-limited quantities and 
    // performing a local half time-step update.

    getBoundsSLOPE(boundL, boundR, slopeCells, u0, w, d);
    localUpdate(boundL, boundR, dx, dt, d);
}


void getBoundsSLOPE(Vector<double>& boundL, Vector<double>& boundR, const Vector<double>& slopeCells,\
                    const Vector<Vector<double> >& u0, double w, const int d){
    
    // Using the 'slopeCells', which are the energy values that surround a cell whose boundary values 
    // we wish to reconstruct, we slope-limit on the conserved variables, extracted from the vector of vector
    // 'u0'. The limiter choice is made in the inputs file.

    double delta, r;
    r = rval(slopeCells[0],slopeCells[1],slopeCells[2]);
    for (int h=0; h<NUM_STATE; h++){
        delta = 0.5*(1+w)*(u0[1][h]-u0[0][h])+0.5*(1-w)*(u0[2][h]-u0[1][h]);
        boundL[h] = u0[1][h]-0.5*getLimiter(r,xsi_R(r))*delta;
        boundR[h] = u0[1][h]+0.5*getLimiter(r,xsi_R(r))*delta;
    }
     
}

void localUpdate(Vector<double>& boundL, Vector<double>& boundR,\
                 double dx, double dt, const int d){
    
    // Perform local half time-step update to achieve second-order in time accuracy.
    
    Vector<double> qL(NUM_STATE);
    Vector<double> qR(NUM_STATE);
    Vector<double> fL(NUM_STATE);
    Vector<double> fR(NUM_STATE);

    qL = {boundL[0],boundL[1],boundL[2],boundL[3],boundL[4]};
    qR = {boundR[0],boundR[1],boundR[2],boundR[3],boundR[4]};
    fL = getEulerFlux(qL,d);
    fR = getEulerFlux(qR,d);
    for (int h=0; h<NUM_STATE; h++){
        boundL[h] = boundL[h]-0.5*(dt/dx)*(fR[h]-fL[h]);  //half time-step information
        boundR[h] = boundR[h]-0.5*(dt/dx)*(fR[h]-fL[h]);  //half time-step information
    }
}

double getLimiter(double r, double rR){

    // This function computes the amount to scale the boundary reconstruction based
    // on the choice of limiter made, and the slope ratio calculated via the energy values
    // bounding the cell we wish to reconstruct.

    if (enLimiter == 0) {     //Superbee

        if (r<=0){
            return 0;
        }
        else if ((r>0)&&(r<=0.5)){
            return 2*r;
        }
        else if ((r>0.5)&&(r<=1)){
            return 1;
        }
        else {
            if ((r<rR)&&(r<2)){
                return r;
            }
            else if ((rR<=r)&&(rR<2)){
                return rR;
            }
            else{
                return 2;
            }
        }

    }

    else if (enLimiter == 1){ //vanLeer

        if (r<=0){
            return 0;
        }
        else {
            if ((2*r/(1+r)) < rR){
                return (2*r/(1+r));
            }
            else {
                return rR;
            }
        }

    }

    else if (enLimiter == 2){ //vanAlbada

        if (r<=0){
            return 0;
        }
        else {
            if (((1+r)*r/(1+r*r))<rR){
                return ((1+r)*r/(1+r*r));
            }
            else {
                return rR;
            }
        }
    }

    else {                  //Minbee

        if (r<=0){
            return 0;
        }
        else if ((r>0)&&(r<=1)){
            return r;
        }
        else {
            if (1<rR){
                return 1;
            }
            else {
                return rR;
            }
        }

    }

}

double rval(double left, double center, double right){
    double num, denom;
    num = center - left;
    denom = right - center;
    
    // In case a division involving zero occurs, we assign arbitrarily 
    // small values to the numerator and denominator.

    if (fabs(num) < 1e-12){
        if (num < 0){
            num = -1e-12;
        }
        else {
            num = 1e-12;
        }
    }
    if (fabs(denom) < 1e-12){
        if (denom < 0){
            denom = -1e-12;
        }
        else {
            denom = 1e-12;
        }
    }
    return num/denom;
}

double xsi_R(double r){
    double val = 2/(1+r);
    return val;
}
