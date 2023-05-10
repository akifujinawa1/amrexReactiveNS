#include "eulerFunc.H"
#include "diffusionFunc.H"
#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>

#include <cmath>
#include <iostream>
#include <algorithm> 

extern int enIC;
extern double Gamma;
extern int NUM_STATE;
extern const int spacedim;
extern double R;
extern double M_O2;
extern double M_N2;

using namespace amrex;

// Use this file to write functions required for thermodynamic and transport calculations

double hO2(const double& Tg){

    double a1,a2,a3,a4,a5,a6,a7,b1,h;
        
    if (Tg <= 1000)
    {
        a1 = -3.425563420*1e4;
        a2 = 4.847000970*1e2;
        a3 = 1.119010961;
        a4 = 4.293889240*1e-3;
        a5 = -6.836300520*1e-7;
        a6 = -2.023372700*1e-9;
        a7 = 1.039040018*1e-12;
        b1 = -3.391454870*1e3;
    }
    else if (Tg <= 6000)
    {
        a1 = -1.037939022e+06;
        a2 = 2.344830282e+03;
        a3 = 1.819732036e+00;
        a4 = 1.267847582e-03;
        a5 = -2.188067988e-07;
        a6 = 2.053719572e-11;
        a7 = -8.193467050e-16;
        b1 = -1.689010929e+04;
    }
    else 
    {
        a1 = 4.975294300e+08;
        a2 = -2.866106874e+05;
        a3 = 6.690352250e+01;
        a4 = -6.169959020e-03;
        a5 = 3.016396027e-07;
        a6 = -7.421416600e-12;
        a7 = 7.278175770e-17;
        b1 = 2.293554027e+06;
    }

    h = (-a1/pow(Tg,2.0) + a2*log(Tg)/Tg + a3 + a4*Tg/2.0 + a5*pow(Tg,2.0)/3.0 + a6*pow(Tg,3.0)/4.0 + a7*pow(Tg,4.0)/5.0 + b1/Tg)*R*Tg;
    return h;

}

double hN2(const double& Tg){

    double a1,a2,a3,a4,a5,a6,a7,b1,h;

    if (Tg <= 1000)
    {
        a1 = 2.210371497e+04;
        a2 = -3.818461820e+02;
        a3 = 6.082738360e+00;
        a4 = -8.530914410e-03;
        a5 = 1.384646189e-05;
        a6 = -9.625793620e-09;
        a7 = 2.519705809e-12;
        b1 = 7.108460860e+02;
    }
    else if (Tg <= 6000)
    {
        a1 = 5.877124060e+05;
        a2 = -2.239249073e+03;
        a3 = 6.066949220e+00;
        a4 = -6.139685500e-04;
        a5 = 1.491806679e-07;
        a6 = -1.923105485e-11;
        a7 = 1.061954386e-15;
        b1 = 1.283210415e+04;
    }
    else 
    {
        a1 = 8.310139160e+08;
        a2 = -6.420733540e+05;
        a3 = 2.020264635e+02;
        a4 = -3.065092046e-02;
        a5 = 2.486903333e-06;
        a6 = -9.705954110e-11;
        a7 = 1.437538881e-15;
        b1 = 4.938707040e+06;
    }

    h = (-a1/pow(Tg,2.0) + a2*log(Tg)/Tg + a3 + a4*Tg/2.0 + a5*pow(Tg,2.0)/3.0 + a6*pow(Tg,3.0)/4.0 + a7*pow(Tg,4.0)/5.0 + b1/Tg)*R*Tg;
    return h;

}

// double hO2(const double& Tg){

//     // Returns sensible enthalpy of O2 gas as a function of gas-phase temperature
//     // in units J/mol, based on the Shomate equation

//     double a,b,c,d,e,f,h,hval;

//     if (Tg < 700){
//         a = 31.32234;
//         b = -20.23531;
//         c = 57.86644;
//         d = -36.50624;
//         e = -0.007374;
//         f = -8.903471;
//         h = 0.0;
//     }
//     else if (Tg < 2000){
//         a = 30.03235;
//         b = 8.772972;
//         c = -3.988133;
//         d = 0.788313;
//         e = -0.741599;
//         f = -11.32468;
//         h = 0.0;
//     }
//     else {
//         a = 20.91111;
//         b = 10.72071;
//         c = -2.020498;
//         d = 0.146449;
//         e = 9.245722;
//         f = 5.337651;
//         h = 0.0;
//     }

//     double Tref = Tg/1000;
//     hval = 1e3*(a*Tref + b*Tref*Tref*0.5 + c*Tref*Tref*Tref*(1.0/3.0) + d*Tref*Tref*Tref*Tref*0.25 - e/Tref + f - h);
//     return hval;
// }

// double hN2(const double& Tg){

//     // Returns sensible enthalpy of O2 gas as a function of gas-phase temperature
//     // in units J/mol, based on the Shomate equation

//     double a,b,c,d,e,f,h,hval;

//     if (Tg < 500){
//         a = 28.98641;
//         b = 1.853978;
//         c = -9.647459;
//         d = 16.63537;
//         e = 0.000117;
//         f = -8.671914;
//         h = 0.0;
//     }
//     else if (Tg < 2000){
//         a = 19.50583;
//         b = 19.88705;
//         c = -8.598535;
//         d = 1.369784;
//         e = 0.527601;
//         f = -4.935202;
//         h = 0.0;
//     }
//     else {
//         a = 35.51872;
//         b = 1.128728;
//         c = -0.196103;
//         d = 0.014662;
//         e = -18.97091;
//         f = 224.9810;
//         h = 0.0;
//     }

//     double Tref = Tg/1000;
//     hval = 1e3*(a*Tref + b*Tref*Tref*0.5 + c*Tref*Tref*Tref*(1.0/3.0) + d*Tref*Tref*Tref*Tref*0.25 - e/Tref + f - h);
//     return hval;
// }

double cpO2(const double& Tg){

    double a1,a2,a3,a4,a5,a6,a7,b1,cp;
        
    if (Tg < 1000)
    {
        a1 = -3.425563420*1e4;
        a2 = 4.847000970*1e2;
        a3 = 1.119010961;
        a4 = 4.293889240*1e-3;
        a5 = -6.836300520*1e-7;
        a6 = -2.023372700*1e-9;
        a7 = 1.039040018*1e-12;
        b1 = -3.391454870*1e3;
    }
    else if (Tg < 6000)
    {
        a1 = -1.037939022e+06;
        a2 = 2.344830282e+03;
        a3 = 1.819732036e+00;
        a4 = 1.267847582e-03;
        a5 = -2.188067988e-07;
        a6 = 2.053719572e-11;
        a7 = -8.193467050e-16;
        b1 = -1.689010929e+04;
    }
    else 
    {
        a1 = 4.975294300e+08;
        a2 = -2.866106874e+05;
        a3 = 6.690352250e+01;
        a4 = -6.169959020e-03;
        a5 = 3.016396027e-07;
        a6 = -7.421416600e-12;
        a7 = 7.278175770e-17;
        b1 = 2.293554027e+06;
    }

    cp = (a1*pow(Tg,-2) + a2/Tg + a3 + a4*Tg + a5*pow(Tg,2.0)  + a6*pow(Tg,3.0) + a7*pow(Tg,4.0))*R;
    return cp;

}

double cpN2(const double& Tg){

    double a1,a2,a3,a4,a5,a6,a7,b1,cp;

    if (Tg < 1000)
    {
        a1 = 2.210371497e+04;
        a2 = -3.818461820e+02;
        a3 = 6.082738360e+00;
        a4 = -8.530914410e-03;
        a5 = 1.384646189e-05;
        a6 = -9.625793620e-09;
        a7 = 2.519705809e-12;
        b1 = 7.108460860e+02;
    }
    else if (Tg < 6000)
    {
        a1 = 5.877124060e+05;
        a2 = -2.239249073e+03;
        a3 = 6.066949220e+00;
        a4 = -6.139685500e-04;
        a5 = 1.491806679e-07;
        a6 = -1.923105485e-11;
        a7 = 1.061954386e-15;
        b1 = 1.283210415e+04;
    }
    else 
    {
        a1 = 8.310139160e+08;
        a2 = -6.420733540e+05;
        a3 = 2.020264635e+02;
        a4 = -3.065092046e-02;
        a5 = 2.486903333e-06;
        a6 = -9.705954110e-11;
        a7 = 1.437538881e-15;
        b1 = 4.938707040e+06;
    }

    cp = (a1*pow(Tg,-2) + a2/Tg + a3 + a4*Tg + a5*pow(Tg,2.0)  + a6*pow(Tg,3.0) + a7*pow(Tg,4.0))*R;
    return cp;

}

double cpMix(const double& cpO2, const double& cpN2, const double& Y_O2, const double& Y_N2){
    double cp = Y_O2*cpO2/M_O2 + Y_N2*cpN2/M_N2;
    return cp;
}

double kO2(const double& Tg){
    double A1,A2,B1,B2,C1,C2,D1,D2,k;
    A1 = 0.77238828; A2 = 0.90875998;
    B1 = 6.9293259;  B2 = 289.86028;
    C1 = -5900.8518; C2 = -79180.433;
    D1 = 1.2202965;  D2 = 0.068622859;
    if ((Tg > 200) && (Tg <= 1000))
    {
        k = 1.0e-4*exp(A1*log(Tg) + B1/Tg + C1/pow(Tg,2.0) + D1);
    }
    else 
    {
        k = 1.0e-4*exp(A2*log(Tg) + B2/Tg + C2/pow(Tg,2.0) + D2);
    }
    return k;
}

double kN2(const double& Tg){
    double A1,A2,B1,B2,C1,C2,D1,D2,k;
    A1 = 0.85372829; A2 = 0.88506520;
    B1 = 105.18665;  B2 = 134.69656;
    C1 = -12299.753; C2 = -11386.420;
    D1 = 0.48299104; D2 = 0.23610008;
    if ((Tg > 200) && (Tg <= 1000))
    {
        k = 1.0e-4*exp(A1*log(Tg) + B1/Tg + C1/pow(Tg,2.0) + D1);
    }
    else 
    {
        k = 1.0e-4*exp(A2*log(Tg) + B2/Tg + C2/pow(Tg,2.0) + D2);
    }
    return k;
}

double kMix(const double& kO2, const double& kN2, const double& Y_O2, const double& Y_N2){

    double X_O2,X_N2,k;

    X_O2 = (Y_O2/M_O2)/(Y_O2/M_O2+Y_N2/M_N2);
    X_N2 = (Y_N2/M_N2)/(Y_O2/M_O2+Y_N2/M_N2);
    k = 0.5*((X_O2*kO2+X_N2*kN2) + 1/(X_O2/kO2+X_N2/kN2));

    return k;
}

double muO2(const double& Tg){
    double A1,A2,B1,B2,C1,C2,D1,D2,mu;
    A1 = 0.60916180; A2 = 0.72216486;
    B1 = -52.244847;  B2 = 175.50839;
    C1 = -599.74009; C2 = -57974.816;
    D1 = 2.0410801; D2 = 1.0901044;
    if ((Tg > 200) && (Tg <= 1000))
    {
        mu = 1.0e-7*exp(A1*log(Tg) + B1/Tg + C1/pow(Tg,2.0) + D1);
    }
    else 
    {
        mu = 1.0e-7*exp(A2*log(Tg) + B2/Tg + C2/pow(Tg,2.0) + D2);
    }
    return mu;
}

double muN2(const double& Tg){
    double A1,A2,B1,B2,C1,C2,D1,D2,mu;
    A1 = 0.62526577; A2 = 0.87395209;
    B1 = -31.779652;  B2 = 561.52222;
    C1 = -1640.7983; C2 = 173948.09;
    D1 = 1.7454992;  D2 = -0.39335958;
    if ((Tg > 200) && (Tg <= 1000))
    {
        mu = 1.0e-7*exp(A1*log(Tg) + B1/Tg + C1/pow(Tg,2.0) + D1);
    }
    else 
    {
        mu = 1.0e-7*exp(A2*log(Tg) + B2/Tg + C2/pow(Tg,2.0) + D2);
    }
    return mu;
}


double Hgas(const double& rho, const double& YO2, const double& YN2, const double& Tg){
    double H = rho*(YO2*hO2(Tg)/M_O2 + YN2*hN2(Tg)/M_N2);
    return H;
}


double getMavg(const double& YO2, const double& YN2){
    double Mavg = 1.0/(YO2/M_O2+YN2/M_N2);
    return Mavg;
}


double Tg(const double& rho, const double& u, const double& v, \
          const double& YO2, const double& YN2, const double& ener){
    double Tg0=300,tol=1e-6,error=1;
    double dE,E,Eest,vT,Mavg,cp,Tgn,iter=0;

    Mavg = getMavg(YO2,YN2);
    vT = sqrt(u*u+v*v);
    E = ener - 0.5*rho*vT*vT;    // volumetric enthalpy of gas

    while (error > tol)
    {
        double enthalpy = Hgas(rho,YO2,YN2,Tg0);
        double pressure = rho*(R/Mavg)*Tg0;
        dE    = E - (enthalpy-pressure);
        // dE    = E - (Hgas(rho,YO2,YN2,Tg0)-rho*(R/Mavg)*Tg0); //
        cp    = rho*(YO2*cpO2(Tg0)/M_O2+YN2*cpN2(Tg0)/M_N2);
        Tgn   = Tg0 + dE/cp;
        error = std::fabs(Tgn-Tg0);
        // std::cout << "ener: " << ener << ", enthalpy guess: " << enthalpy << \
        // ", p guess: " << pressure << ", Tgn: " << Tgn << ", error: " << error << std::endl;
        Tg0 = Tgn;
        iter += 1;
        if (iter > 20){
            Abort("Gas temperature iteration does not converge"); 
        }
    }

    return Tgn;

}

Vector<double> getMoleFractions(const double& YO2, const double& YN2){
    Vector<double> xGas(2);
    xGas[0] = (YO2/M_O2)/(YO2/M_O2 + YN2/M_N2);
    xGas[1] = (YN2/M_N2)/(YO2/M_O2 + YN2/M_N2);
    return xGas;    
}

Vector<double> getMassFractions(const double& XO2, const double& XN2){
    Vector<double> yGas(2);
    yGas[0] = (XO2*M_O2)/(XO2*M_O2 + XN2*M_N2);
    yGas[1] = (XN2*M_N2)/(XO2*M_O2 + XN2*M_N2);
    return yGas;    
}

double muMix(const double& muO2, const double& muN2, const double& YO2, const double& YN2){
    double phiO2N2,phiN2O2,XO2,XN2,mu;
    Vector<double> xGas(2);
    xGas = getMoleFractions(YO2,YN2);
    XO2 = xGas[0];
    XN2 = xGas[1];
    phiO2N2 = pow((1/sqrt(8.0))*sqrt(1+(M_O2/M_N2))*(1+sqrt(muO2/muN2)*pow(M_N2/M_O2,0.25)),2);
    phiN2O2 = pow((1/sqrt(8.0))*sqrt(1+(M_N2/M_O2))*(1+sqrt(muN2/muO2)*pow(M_O2/M_N2,0.25)),2);
    mu = XO2*muO2/(XO2+XN2*phiO2N2) + XN2*muN2/(XN2+XO2*phiN2O2);
    return mu;
}


