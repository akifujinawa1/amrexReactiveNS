// This file contains various functions that are used to compute the source terms
// for the augmented Euler equations.

#include "source.H"
#include "eulerFunc.H"
#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

#include <math.h>
#include <iostream>
#include <algorithm> 

extern int source;
extern int NUM_STATE;

extern double Gamma;
extern double R;
extern double M;
extern double A;
extern double T0;
extern double Ea;
extern double q;

using namespace amrex;

void updateSource(MultiFab& Sborder, Vector<double>& q, Vector<double>& sourceVec, const double &dt, const int &source){

    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Dim3 lo = lbound(bx);
        const Dim3 hi = ubound(bx);
        const auto& arr = Sborder.array(mfi);

        for(int k = lo.z; k <= hi.z; k++)
        {
            for(int j = lo.y; j <= hi.y; j++)
            {
                for(int i = lo.x; i <= hi.x; i++)
                {
                    for(int h = 0; h < NUM_STATE; h++)
                    {
                        q[h] = arr(i,j,k,h);
                    }
                    getSource(q,sourceVec,source,dt);
                    for(int h = 0; h < NUM_STATE; h++)
                    {
                        arr(i,j,k,h) += sourceVec[h];
                    }
                    // apply limit on lambda here to avoid unboundedness
                    double lambda = arr(i,j,k,4)/arr(i,j,k,0);
                    if (lambda < 1e-5)
                    {
                        arr(i,j,k,4) = 0.0;
                    }
                }
            }
        }
    }
}

void getSource(const Vector<double>& q, Vector<double>& sourceVec, const int& source, const double& dt){
    
    // We update a vector (passed by ref) with the source term contributions.

    // The choice of ODE solver is made in the inputs file. If source == 1, explicit Euler, if 
    // source == 2, RK4.

    switch (source)
    {
        case 1: //explicit euler
        {
            EulerEx(sourceVec,q,dt);
            break;
        }
        case 2: //rk4
        {
            rk4(sourceVec,q,dt);
            break;
        }
    }
}

void EulerEx(Vector<double>& source, const Vector<double>& q, const double& dt){

    // This is the explicit Euler method to solve the ODE. Mainly implemented
    // for code validation, use RK4 instead for better temporal accuracy.

    double rho, u, ener, p, lambda, T;

    rho    = q[0];
    u      = q[1]/rho;
    ener   = q[3];
    lambda = q[4]/rho;
    p      = pressure(rho,u,0,ener);
    T      = p/(rho*(R/M));

    source[0]=0;
    source[1]=0;
    source[2]=0;
    source[3]=dt*rho*getEnergySource(rho,lambda,T);     
    source[4]=dt*getReactiveSource(rho,lambda,T); // conserved variable is rhoY, source term ODE is dY/dt

}

void rk4(Vector<double>& source, const Vector<double>& q, const double& dt){

    // This is the fourth-order Runge-Kutta method to solve the ODE. While exact source of the
    // algorithm is unknown, history of Runge-Kutta methods can be found in:
    // J. C. Butcher, A history of runge-kutta methods, Applied numerical mathematics 20 (3) (1996) 247–260
    // This method is preferred over the explicit Euler method for its superior temporal accuracy.


    double K1e,K2e,K3e,K4e,K1r,K2r,K3r,K4r,newE,newP,newT;
    double rho, u, ener, p, lambda, T, newLambda;
        
    rho    = q[0];
    u      = q[1]/rho;
    ener   = q[3];
    lambda = q[4]/rho;
    p      = pressure(rho,u,0,ener);
    T      = p/(rho*(R/M));

    K1e = dt*getEnergySource(rho,lambda,T);
    K1r = dt*getReactiveSource(rho,lambda,T);
    newLambda = lambda+0.5*K1r;
    newE = ener+0.5*K1e;
    newP = pressure(rho,u,0,newE);
    newT = newP/(rho*(R/M));

    K2e = dt*getEnergySource(rho,newLambda,newT);
    K2r = dt*getReactiveSource(rho,newLambda,newT);
    newLambda = lambda+0.5*K2r;
    newE = ener+0.5*K2e;
    newP = pressure(rho,u,0,newE);
    newT = newP/(rho*(R/M));

    K3e = dt*getEnergySource(rho,newLambda,newT);
    K3r = dt*getReactiveSource(rho,newLambda,newT);
    newLambda = lambda+K3r;
    newE = ener+K3e;
    newP = pressure(rho,u,0,newE);
    newT = newP/(rho*(R/M));

    K4e = dt*getEnergySource(rho,newLambda,newT);
    K4r = dt*getReactiveSource(rho,newLambda,newT);

    source[0]=0;
    source[1]=0;
    source[2]=0;
    source[3]=(1.0/6.0)*(K1e+2.0*K2e+2.0*K3e+K4e);
    source[4]=(1.0/6.0)*rho*(K1r+2.0*K2r+2.0*K3r+K4r);
}


double getEnergySource(const double& rho, const double& lambda, const double& T){

    // Computes energy source term based on the source term in the formulation found in:
    // Michael Liberman, Cheng Wang, Chengeng Qian & JianNan Liu (2019)
    // Influence of chemical kinetics on spontaneous waves and detonation initiation in highly
    // reactive and low reactive mixtures, Combustion Theory and Modelling, 23:3, 467-495, DOI:
    // 10.1080/13647830.2018.1551578

    double source;
    source = rho*(q*R*T0/M)*A*rho*lambda*exp(-Ea*T0/T); // rho added for dimensional agreement (source = J/m^3/s)
    return source;
}

double getReactiveSource(const double& rho, const double& lambda, const double& T){

    // Computes reactive source term based on the source term in the formulation found in:
    // Michael Liberman, Cheng Wang, Chengeng Qian & JianNan Liu (2019)
    // Influence of chemical kinetics on spontaneous waves and detonation initiation in highly
    // reactive and low reactive mixtures, Combustion Theory and Modelling, 23:3, 467-495, DOI:
    // 10.1080/13647830.2018.1551578

    double source;
    source = -A*rho*lambda*exp(-Ea*T0/T);
    return source;
}