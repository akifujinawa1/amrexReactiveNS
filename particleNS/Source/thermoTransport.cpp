#include "eulerFunc.H"
#include "diffusionFunc.H"
#include "AmrLevelAdv.H"
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
extern double M_O2, M_N2, M_Fe, M_FeO;
extern double kb;
extern double sigmaO2N2, sigmaO2Fe, sigmaO2FeO, sigmaN2Fe, sigmaN2FeO, sigmaFeFeO;
extern double epsO2N2, epsO2Fe, epsO2FeO, epsN2Fe, epsN2FeO, epsFeFeO;
extern double M_O2N2, M_O2Fe, M_O2FeO, M_N2Fe, M_N2FeO, M_FeFeO;

using namespace amrex;

// Use this file to write functions required for thermodynamic and transport calculations

// double hO2(const double& Tg){

//     double a1,a2,a3,a4,a5,a6,a7,b1,h;
        
//     if (Tg <= 1000)
//     {
//         a1 = -3.425563420*1e4;
//         a2 = 4.847000970*1e2;
//         a3 = 1.119010961;
//         a4 = 4.293889240*1e-3;
//         a5 = -6.836300520*1e-7;
//         a6 = -2.023372700*1e-9;
//         a7 = 1.039040018*1e-12;
//         b1 = -3.391454870*1e3;
//     }
//     else if (Tg <= 6000)
//     {
//         a1 = -1.037939022e+06;
//         a2 = 2.344830282e+03;
//         a3 = 1.819732036e+00;
//         a4 = 1.267847582e-03;
//         a5 = -2.188067988e-07;
//         a6 = 2.053719572e-11;
//         a7 = -8.193467050e-16;
//         b1 = -1.689010929e+04;
//     }
//     else 
//     {
//         a1 = 4.975294300e+08;
//         a2 = -2.866106874e+05;
//         a3 = 6.690352250e+01;
//         a4 = -6.169959020e-03;
//         a5 = 3.016396027e-07;
//         a6 = -7.421416600e-12;
//         a7 = 7.278175770e-17;
//         b1 = 2.293554027e+06;
//     }

//     h = (-a1/pow(Tg,2.0) + a2*log(Tg)/Tg + a3 + a4*Tg/2.0 + a5*pow(Tg,2.0)/3.0 + a6*pow(Tg,3.0)/4.0 + a7*pow(Tg,4.0)/5.0 + b1/Tg)*R*Tg;
//     return h;

// }

// double hN2(const double& Tg){

//     double a1,a2,a3,a4,a5,a6,a7,b1,h;

//     if (Tg <= 1000)
//     {
//         a1 = 2.210371497e+04;
//         a2 = -3.818461820e+02;
//         a3 = 6.082738360e+00;
//         a4 = -8.530914410e-03;
//         a5 = 1.384646189e-05;
//         a6 = -9.625793620e-09;
//         a7 = 2.519705809e-12;
//         b1 = 7.108460860e+02;
//     }
//     else if (Tg <= 6000)
//     {
//         a1 = 5.877124060e+05;
//         a2 = -2.239249073e+03;
//         a3 = 6.066949220e+00;
//         a4 = -6.139685500e-04;
//         a5 = 1.491806679e-07;
//         a6 = -1.923105485e-11;
//         a7 = 1.061954386e-15;
//         b1 = 1.283210415e+04;
//     }
//     else 
//     {
//         a1 = 8.310139160e+08;
//         a2 = -6.420733540e+05;
//         a3 = 2.020264635e+02;
//         a4 = -3.065092046e-02;
//         a5 = 2.486903333e-06;
//         a6 = -9.705954110e-11;
//         a7 = 1.437538881e-15;
//         b1 = 4.938707040e+06;
//     }

//     h = (-a1/pow(Tg,2.0) + a2*log(Tg)/Tg + a3 + a4*Tg/2.0 + a5*pow(Tg,2.0)/3.0 + a6*pow(Tg,3.0)/4.0 + a7*pow(Tg,4.0)/5.0 + b1/Tg)*R*Tg;
//     return h;

// }

double hO2(const double& Tg){

    // Returns sensible enthalpy of O2 gas as a function of gas-phase temperature
    // in units J/mol, based on the Shomate equation

    double a,b,c,d,e,f,h,hval;

    if (Tg < 700){
        a = 31.32234;
        b = -20.23531;
        c = 57.86644;
        d = -36.50624;
        e = -0.007374;
        f = -8.903471;
        h = 0.0;
    }
    else if (Tg < 2000){
        a = 30.03235;
        b = 8.772972;
        c = -3.988133;
        d = 0.788313;
        e = -0.741599;
        f = -11.32468;
        h = 0.0;
    }
    else {
        a = 20.91111;
        b = 10.72071;
        c = -2.020498;
        d = 0.146449;
        e = 9.245722;
        f = 5.337651;
        h = 0.0;
    }

    double Tref = Tg/1000;
    hval = 1e3*(a*Tref + b*Tref*Tref*0.5 + c*Tref*Tref*Tref*(1.0/3.0) + d*Tref*Tref*Tref*Tref*0.25 - e/Tref + f - h);
    return hval;
}

double hN2(const double& Tg){

    // Returns sensible enthalpy of O2 gas as a function of gas-phase temperature
    // in units J/mol, based on the Shomate equation

    double a,b,c,d,e,f,h,hval;

    if (Tg < 500){
        a = 28.98641;
        b = 1.853978;
        c = -9.647459;
        d = 16.63537;
        e = 0.000117;
        f = -8.671914;
        h = 0.0;
    }
    else if (Tg < 2000){
        a = 19.50583;
        b = 19.88705;
        c = -8.598535;
        d = 1.369784;
        e = 0.527601;
        f = -4.935202;
        h = 0.0;
    }
    else {
        a = 35.51872;
        b = 1.128728;
        c = -0.196103;
        d = 0.014662;
        e = -18.97091;
        f = 224.9810;
        h = 0.0;
    }

    double Tref = Tg/1000;
    hval = 1e3*(a*Tref + b*Tref*Tref*0.5 + c*Tref*Tref*Tref*(1.0/3.0) + d*Tref*Tref*Tref*Tref*0.25 - e/Tref + f - h);
    return hval;
}

double cpO2(const double& Tg){

    double a1,a2,a3,a4,a5,a6,a7,b1,cp;
        
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

    cp = (a1*pow(Tg,-2) + a2/Tg + a3 + a4*Tg + a5*pow(Tg,2.0)  + a6*pow(Tg,3.0) + a7*pow(Tg,4.0))*R;
    return cp;

}

double cpN2(const double& Tg){

    double a1,a2,a3,a4,a5,a6,a7,b1,cp;

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

    cp = (a1*pow(Tg,-2.0) + a2/Tg + a3 + a4*Tg + a5*pow(Tg,2.0)  + a6*pow(Tg,3.0) + a7*pow(Tg,4.0))*R;
    return cp;

}

double cpMix(const double& cpO2, const double& cpN2, const double& Y_O2, const double& Y_N2){
    double cp = Y_O2*cpO2/M_O2 + Y_N2*cpN2/M_N2;
    return cp;
}

double cpMixFe(const double& cpO2, const double& cpN2, const double& cpFe, const double& cpFeO, \
               const double& Y_O2, const double& Y_N2, const double& Y_Fe, const double& Y_FeO){
    double cp = Y_O2*cpO2/M_O2 + Y_N2*cpN2/M_N2 + Y_Fe*cpFe/M_Fe + Y_FeO*cpFeO/M_FeO;
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

double kFe(const double& Tg){
    // note that this function is used for FeO as well, as the same LJ parameters (well-depth and diameter) 
    // are used for gas-phase Fe and FeO. 
    double A1,A2,B1,B2,C1,C2,D1,D2,E1,k;
    A1 = -1.158e-14; A2 = 1.657;
    B1 = 3.539e-11;  B2 = 915.1;
    C1 = -3.975e-8; C2 = -2.502e4;
    D1 = 2.523e-5; D2 = -8.032;
    E1 = -0.001842;
    if ((Tg > 200) && (Tg <= 1000))
    {
        k = A1*pow(Tg,4.0) + B1*pow(Tg,3.0) + C1*pow(Tg,2.0) + D1*Tg + E1;
    }
    else 
    {
        k = 1.0e-4*exp(A2*log(Tg) + B2/Tg + C2/pow(Tg,2.0) + D2);
    }
    return k;
}

double kFeO(const double& Tg){
    // note that this function is used for FeO as well, as the same LJ parameters (well-depth and diameter) 
    // are used for gas-phase Fe and FeO. 
    double A1,A2,B1,B2,C1,C2,D1,D2,E1,k;
    A1 = -1.158e-14; A2 = 1.657;
    B1 = 3.539e-11;  B2 = 915.1;
    C1 = -3.975e-8; C2 = -2.502e4;
    D1 = 2.523e-5; D2 = -8.032;
    E1 = -0.001842;
    if ((Tg > 200) && (Tg <= 1000))
    {
        k = A1*pow(Tg,4.0) + B1*pow(Tg,3.0) + C1*pow(Tg,2.0) + D1*Tg + E1;
    }
    else 
    {
        k = 1.0e-4*exp(A2*log(Tg) + B2/Tg + C2/pow(Tg,2.0) + D2);
    }
    return k;
}

Vector<double> getMoleFractions(const double& Y_O2, const double& Y_N2, const double& Y_Fe, const double& Y_FeO){
    Vector<double> xGas(4);
    xGas[0] = (Y_O2/M_O2)/(Y_O2/M_O2+Y_N2/M_N2+Y_Fe/M_Fe+Y_FeO/M_FeO);
    xGas[1] = (Y_N2/M_N2)/(Y_O2/M_O2+Y_N2/M_N2+Y_Fe/M_Fe+Y_FeO/M_FeO);
    xGas[2] = (Y_Fe/M_Fe)/(Y_O2/M_O2+Y_N2/M_N2+Y_Fe/M_Fe+Y_FeO/M_FeO);
    xGas[3] = (Y_FeO/M_FeO)/(Y_O2/M_O2+Y_N2/M_N2+Y_Fe/M_Fe+Y_FeO/M_FeO);
    return xGas;    
}

Vector<double> getMassFractions(const double& X_O2, const double& X_N2, const double& X_Fe, const double& X_FeO){
    Vector<double> yGas(4);
    yGas[0] = (X_O2*M_O2)/(X_O2*M_O2 + X_N2*M_N2 + X_Fe*M_Fe + X_FeO*M_FeO);
    yGas[1] = (X_N2*M_N2)/(X_O2*M_O2 + X_N2*M_N2 + X_Fe*M_Fe + X_FeO*M_FeO);
    yGas[2] = (X_Fe*M_Fe)/(X_O2*M_O2 + X_N2*M_N2 + X_Fe*M_Fe + X_FeO*M_FeO);
    yGas[3] = (X_FeO*M_FeO)/(X_O2*M_O2 + X_N2*M_N2 + X_Fe*M_Fe + X_FeO*M_FeO);
    return yGas;    
}

double kMix(const double& kO2, const double& kN2, const double& kFe, const double& kFeO,\
            const double& Y_O2, const double& Y_N2, const double& Y_Fe, const double& Y_FeO){

    // Wassiljewa equation with Mason and Saxena mixture rule (similar to viscosity mixture rule)
    // We use this equation for mixture thermal conductivity, since there is a large difference in 
    // molecular weights between the gas species we consider (M_FeO/M_N2 > 2.5)

    // double XO2,XN2,XFe,XFeO,k,a,b,c,d;
    // double phiO2N2,phiO2Fe,phiO2FeO,phiN2O2,phiN2Fe,phiN2FeO;
    // double phiFeO2,phiFeN2,phiFeFeO,phiFeOO2,phiFeON2,phiFeOFe;
    // Vector<double> xGas(gases::ncomps,0);
    // xGas = getMoleFractions(Y_O2,Y_N2,Y_Fe,Y_FeO);
    // XO2 = xGas[0];
    // XN2 = xGas[1];
    // XFe = xGas[2];
    // XFeO= xGas[3];

    // phiO2N2 = (1/sqrt(8.0))*(1/sqrt(1+(M_O2/M_N2)))*pow(1+sqrt(muO2/muN2)*pow(M_N2/M_O2,0.25),2);
    // phiO2Fe = (1/sqrt(8.0))*(1/sqrt(1+(M_O2/M_Fe)))*pow(1+sqrt(muO2/muFe)*pow(M_Fe/M_O2,0.25),2);
    // phiO2FeO= (1/sqrt(8.0))*(1/sqrt(1+(M_O2/M_FeO)))*pow(1+sqrt(muO2/muFeO)*pow(M_FeO/M_O2,0.25),2);

    // phiN2O2 = (muN2/muO2)*(M_O2/M_N2)*phiO2N2;
    // phiN2Fe = (1/sqrt(8.0))*(1/sqrt(1+(M_N2/M_Fe)))*pow(1+sqrt(muN2/muFe)*pow(M_Fe/M_N2,0.25),2);
    // phiN2FeO = (1/sqrt(8.0))*(1/sqrt(1+(M_N2/M_FeO)))*pow(1+sqrt(muN2/muFeO)*pow(M_FeO/M_N2,0.25),2);

    // phiFeO2 = (muFe/muO2)*(M_O2/M_Fe)*phiO2Fe;
    // phiFeN2 = (muFe/muN2)*(M_N2/M_Fe)*phiN2Fe;
    // phiFeFeO= (1/sqrt(8.0))*(1/sqrt(1+(M_Fe/M_FeO)))*pow(1+sqrt(muFe/muFeO)*pow(M_FeO/M_Fe,0.25),2);

    // phiFeOO2 = (muFeO/muO2)*(M_O2/M_FeO)*phiO2FeO;
    // phiFeON2 = (muFeO/muN2)*(M_N2/M_FeO)*phiN2FeO;
    // phiFeOFe = (muFeO/muFe)*(M_Fe/M_FeO)*phiFeFeO;

    // a = XO2*kO2/(XO2+XN2*phiO2N2+XFe*phiO2Fe+XFeO*phiO2FeO);
    // b = XN2*kN2/(XN2+XO2*phiN2O2+XFe*phiN2Fe+XFeO*phiN2FeO);
    // c = XFe*kFe/(XFe+XO2*phiFeO2+XN2*phiFeN2+XFeO*phiFeFeO);
    // d = XFeO*kFeO/(XFeO+XO2*phiFeOO2+XN2*phiFeON2+XFe*phiFeOFe);
    // k = a+b+c+d;

    // return k;

    // Wilke mixture rule

    double XO2,XN2,XFe,XFeO,k;
    
    Vector<double> xGas(gases::ncomps,0);
    xGas = getMoleFractions(Y_O2,Y_N2,Y_Fe,Y_FeO);
    XO2 = xGas[0];
    XN2 = xGas[1];
    XFe = xGas[2];
    XFeO= xGas[3];

    k = 0.5*((XO2*kO2+XN2*kN2+XFe*kFe+XFeO*kFeO) + 1.0/(XO2/kO2+XN2/kN2+XFe/kFe+XFeO/kFeO));
    return k;
}

double kMix_O2N2(const double& kO2, const double& kN2, const double& Y_O2, const double& Y_N2){

    // Wilke mixture rule

    double XO2,XN2,k;
    
    XO2 = (Y_O2/M_O2)/(Y_O2/M_O2+Y_N2/M_N2);
    XN2 = (Y_N2/M_N2)/(Y_O2/M_O2+Y_N2/M_N2);

    k = 0.5*((XO2*kO2+XN2*kN2) + 1.0/(XO2/kO2+XN2/kN2));
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
    C1 = -1640.7983; C2 = -173948.09;
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

double muFe(const double& Tg){
    // note that this function is used for FeO as well, as the same LJ parameters (well-depth and diameter) 
    // are used for gas-phase Fe and FeO. 
    double A1,A2,B1,B2,C1,C2,D1,D2,E1,E2,mu;
    A1 = -1.111e-17; A2 = 1.58e-20;
    B1 = 3.269e-14;  B2 = -2.884e-16;
    C1 = -3.514e-11; C2 = 1.662e-12;
    D1 = 2.735e-8;   D2 = 8.758e-9;
    E1 = -1.309e-6;  E2 = 2.355e-6;
    if ((Tg > 200) && (Tg <= 1000))
    {
        mu = A1*pow(Tg,4.0) + B1*pow(Tg,3.0) + C1*pow(Tg,2.0) + D1*Tg + E1;
    }
    else 
    {
        mu = A2*pow(Tg,4.0) + B2*pow(Tg,3.0) + C2*pow(Tg,2.0) + D2*Tg + E2;
    }
    return mu;
}

double muFeO(const double& Tg){
    // note that this function is used for FeO as well, as the same LJ parameters (well-depth and diameter) 
    // are used for gas-phase Fe and FeO. 
    double A1,A2,B1,B2,C1,C2,D1,D2,E1,E2,mu;
    A1 = -1.111e-17; A2 = 1.58e-20;
    B1 = 3.269e-14;  B2 = -2.884e-16;
    C1 = -3.514e-11; C2 = 1.662e-12;
    D1 = 2.735e-8;   D2 = 8.758e-9;
    E1 = -1.309e-6;  E2 = 2.355e-6;
    if ((Tg > 200) && (Tg <= 1000))
    {
        mu = A1*pow(Tg,4.0) + B1*pow(Tg,3.0) + C1*pow(Tg,2.0) + D1*Tg + E1;
    }
    else 
    {
        mu = A2*pow(Tg,4.0) + B2*pow(Tg,3.0) + C2*pow(Tg,2.0) + D2*Tg + E2;
    }
    return mu;
}

double Hgas(const double& rho, const double& YO2, const double& YN2, const double& Tg){
    // returns gas-phase sensible enthalpy (std H of form. of O2, N2 are zero) of the 
    // binary O2-N2 gas mixture as a function of temperature, YO2, and YN2, in units
    // J/m^3 (scaled by density to make volumetric)
    double H = rho*(YO2*hO2(Tg)/M_O2 + YN2*hN2(Tg)/M_N2);
    return H;
}


double getMavg(const double& YO2, const double& YN2){
    double Mavg = 1.0/(YO2/M_O2+YN2/M_N2);
    return Mavg;
}


double Tg(const double& rho, const double& u, const double& v, \
          const double& YO2, const double& YN2, const double& ener){
    double Tg0=300,tol=1e-4,error=1;
    double dE,E,Eest,vT,Mavg,cp,Tgn,iter=0;

    Mavg = getMavg(YO2,YN2);
    vT = sqrt(u*u+v*v);
    E = ener - 0.5*rho*vT*vT;    // volumetric enthalpy of gas

    while (error > tol)
    {
        double enthalpy = Hgas(rho,YO2,YN2,Tg0);
        double pressure = rho*(R/Mavg)*Tg0;
        dE    = E - (enthalpy-pressure);
        cp    = rho*(YO2*cpO2(Tg0)/M_O2+YN2*cpN2(Tg0)/M_N2);
        Tgn   = Tg0 + dE/cp;
        error = std::fabs(Tgn-Tg0)/std::fabs(Tgn);
        Print() << "rho: " << rho << ", ener: " << ener << ", enthalpy guess: " << enthalpy << \
        ", p guess: " << pressure << ", Tgn: " << Tgn << ", error: " << error << std::endl;
        Tg0 = Tgn;
        iter += 1;
        if (iter > 20){
            Abort("Gas temperature iteration does not converge"); 
        }
    }

    return Tgn;

}

double muMix(const double& muO2, const double& muN2, const double& muFe, const double& muFeO, \
             const double& YO2, const double& YN2, const double& YFe, const double& YFeO){
    // Uses the Herning and Zipperer mixture rule to calculate mixture viscosity based on mole fractions
    // and viscosity of individual gas species
    double XO2,XN2,XFe,XFeO,mu,a,b,c,d;
    double phiO2N2,phiO2Fe,phiO2FeO,phiN2O2,phiN2Fe,phiN2FeO;
    double phiFeO2,phiFeN2,phiFeFeO,phiFeOO2,phiFeON2,phiFeOFe;
    Vector<double> xGas(gases::ncomps,0);
    xGas = getMoleFractions(YO2,YN2,YFe,YFeO);
    XO2 = xGas[0];
    XN2 = xGas[1];
    XFe = xGas[2];
    XFeO= xGas[3];

    phiO2N2 = (1.0/sqrt(8.0))*(1.0/sqrt(1.0+(M_O2/M_N2)))*pow(1.0+sqrt(muO2/muN2)*pow(M_N2/M_O2,0.25),2.0);
    phiO2Fe = (1.0/sqrt(8.0))*(1.0/sqrt(1.0+(M_O2/M_Fe)))*pow(1.0+sqrt(muO2/muFe)*pow(M_Fe/M_O2,0.25),2.0);
    phiO2FeO= (1.0/sqrt(8.0))*(1.0/sqrt(1.0+(M_O2/M_FeO)))*pow(1.0+sqrt(muO2/muFeO)*pow(M_FeO/M_O2,0.25),2.0);

    phiN2O2 = (muN2/muO2)*(M_O2/M_N2)*phiO2N2;
    phiN2Fe = (1.0/sqrt(8.0))*(1.0/sqrt(1.0+(M_N2/M_Fe)))*pow(1.0+sqrt(muN2/muFe)*pow(M_Fe/M_N2,0.25),2.0);
    phiN2FeO = (1.0/sqrt(8.0))*(1.0/sqrt(1.0+(M_N2/M_FeO)))*pow(1.0+sqrt(muN2/muFeO)*pow(M_FeO/M_N2,0.25),2.0);

    phiFeO2 = (muFe/muO2)*(M_O2/M_Fe)*phiO2Fe;
    phiFeN2 = (muFe/muN2)*(M_N2/M_Fe)*phiN2Fe;
    phiFeFeO= (1.0/sqrt(8.0))*(1.0/sqrt(1.0+(M_Fe/M_FeO)))*pow(1.0+sqrt(muFe/muFeO)*pow(M_FeO/M_Fe,0.25),2.0);

    phiFeOO2 = (muFeO/muO2)*(M_O2/M_FeO)*phiO2FeO;
    phiFeON2 = (muFeO/muN2)*(M_N2/M_FeO)*phiN2FeO;
    phiFeOFe = (muFeO/muFe)*(M_Fe/M_FeO)*phiFeFeO;

    a = XO2*muO2/(XO2+XN2*phiO2N2+XFe*phiO2Fe+XFeO*phiO2FeO);
    b = XN2*muN2/(XN2+XO2*phiN2O2+XFe*phiN2Fe+XFeO*phiN2FeO);
    c = XFe*muFe/(XFe+XO2*phiFeO2+XN2*phiFeN2+XFeO*phiFeFeO);
    d = XFeO*muFeO/(XFeO+XO2*phiFeOO2+XN2*phiFeON2+XFe*phiFeOFe);
    mu = a+b+c+d;

    // mu = XO2*muO2/(XO2+XN2*phiO2N2) + XN2*muN2/(XN2+XO2*phiN2O2);
    // mu = (muFe*XFe*M_Fe + muFeO*XFeO*M_FeO + muO2*XO2*M_O2 + muN2*XN2*M_N2)/(XFe*M_Fe + XFeO*M_FeO + XO2*M_O2 + XN2*M_N2);
    return mu;
}

double muMix_O2N2(const double& muO2, const double& muN2, const double& YO2, const double& YN2){
    // Uses the Herning and Zipperer mixture rule to calculate mixture viscosity based on mole fractions
    // and viscosity of individual gas species
    double XO2,XN2,mu,a,b;
    double phiO2N2,phiN2O2;
    
    XO2 = (YO2/M_O2)/(YO2/M_O2+YN2/M_N2);
    XN2 = (YN2/M_N2)/(YO2/M_O2+YN2/M_N2);

    phiO2N2 = (1.0/sqrt(8.0))*(1.0/sqrt(1.0+(M_O2/M_N2)))*pow(1.0+sqrt(muO2/muN2)*pow(M_N2/M_O2,0.25),2.0);
    phiN2O2 = (muN2/muO2)*(M_O2/M_N2)*phiO2N2;
    
    a = XO2*muO2/(XO2+XN2*phiO2N2);
    b = XN2*muN2/(XN2+XO2*phiN2O2);
    mu = a+b;

    // mu = XO2*muO2/(XO2+XN2*phiO2N2) + XN2*muN2/(XN2+XO2*phiN2O2);
    // mu = (muFe*XFe*M_Fe + muFeO*XFeO*M_FeO + muO2*XO2*M_O2 + muN2*XN2*M_N2)/(XFe*M_Fe + XFeO*M_FeO + XO2*M_O2 + XN2*M_N2);
    return mu;
}


double omega(const double& T, const double& eps){
    // calculates the collision integral based on method of Neufield et al. (1972)
    double Tstar = T*kb/eps;
    double a,b,c,d,e,f,g,h;
    a = 1.06036;
    b = 0.15610;
    c = 0.19300;
    d = 0.47635;
    e = 1.03587;
    f = 1.52996;
    g = 1.76474;
    h = 3.89411;
    double omega_val = a/pow(Tstar,b) + c/exp(d*Tstar) + e/exp(f*Tstar) + g/exp(h*Tstar);
    // std::cout << "Tstar: " << Tstar << " omega: " << omega_val << std::endl;
    return omega_val;
}

Vector<double> getMixDiffCoeffs(const double& T, const double& p, const double& YO2, \
                                const double& YN2, const double& YFe, const double& YFeO){
    
    double XO2, XN2, XFe, XFeO, p_bar;
    double DO2N2, DO2Fe, DO2FeO, DN2Fe, DN2FeO, DFeFeO, M_mean;
    Vector<double> mixDiffCoeffs(gases::ncomps,0);
    Vector<double> xGas(gases::ncomps,0);

    xGas = getMoleFractions(YO2,YN2,YFe,YFeO);
    XO2 = xGas[0];
    XN2 = xGas[1];
    XFe = xGas[2];
    XFeO= xGas[3];

    M_mean = XO2*M_O2+XN2*M_N2+XFe*M_Fe+XFeO*M_FeO;

    p_bar = p/1.0e5;

    // get binary diffusion coefficients of each binary gas pair:    0.00266

    DO2N2  = 0.00266*1e-4*pow(T,3.0/2.0)/(p_bar*sqrt(M_O2N2)*sigmaO2N2*sigmaO2N2*omega(T,epsO2N2));
    DO2Fe  = 0.00266*1e-4*pow(T,3.0/2.0)/(p_bar*sqrt(M_O2Fe)*sigmaO2Fe*sigmaO2Fe*omega(T,epsO2Fe));
    DO2FeO = 0.00266*1e-4*pow(T,3.0/2.0)/(p_bar*sqrt(M_O2FeO)*sigmaO2FeO*sigmaO2FeO*omega(T,epsO2FeO));
    DN2Fe  = 0.00266*1e-4*pow(T,3.0/2.0)/(p_bar*sqrt(M_N2Fe)*sigmaN2Fe*sigmaN2Fe*omega(T,epsN2Fe));
    DN2FeO = 0.00266*1e-4*pow(T,3.0/2.0)/(p_bar*sqrt(M_N2FeO)*sigmaN2FeO*sigmaN2FeO*omega(T,epsN2FeO));
    DFeFeO = 0.00266*1e-4*pow(T,3.0/2.0)/(p_bar*sqrt(M_FeFeO)*sigmaFeFeO*sigmaFeFeO*omega(T,epsFeFeO));

    // DO2N2  = (3.03-0.98/sqrt(M_O2N2))*1e-7*pow(T,3.0/2.0)/(p_bar*sqrt(M_O2N2)*sigmaO2N2*sigmaO2N2*omega(T,epsO2N2));
    // DO2Fe  = (3.03-0.98/sqrt(M_O2Fe))*1e-7*pow(T,3.0/2.0)/(p_bar*sqrt(M_O2Fe)*sigmaO2Fe*sigmaO2Fe*omega(T,epsO2Fe));
    // DO2FeO = (3.03-0.98/sqrt(M_O2FeO))*1e-7*pow(T,3.0/2.0)/(p_bar*sqrt(M_O2FeO)*sigmaO2FeO*sigmaO2FeO*omega(T,epsO2FeO));
    // DN2Fe  = (3.03-0.98/sqrt(M_N2Fe))*1e-7*pow(T,3.0/2.0)/(p_bar*sqrt(M_N2Fe)*sigmaN2Fe*sigmaN2Fe*omega(T,epsN2Fe));
    // DN2FeO = (3.03-0.98/sqrt(M_N2FeO))*1e-7*pow(T,3.0/2.0)/(p_bar*sqrt(M_N2FeO)*sigmaN2FeO*sigmaN2FeO*omega(T,epsN2FeO));
    // DFeFeO = (3.03-0.98/sqrt(M_FeFeO))*1e-7*pow(T,3.0/2.0)/(p_bar*sqrt(M_FeFeO)*sigmaFeFeO*sigmaFeFeO*omega(T,epsFeFeO));

    mixDiffCoeffs[gases::O2]  = ((M_mean-XO2*M_O2)/M_mean)/(XN2/DO2N2 + XFe/DO2Fe + XFeO/DO2FeO);
    mixDiffCoeffs[gases::N2]  = ((M_mean-XN2*M_N2)/M_mean)/(XO2/DO2N2 + XFe/DN2Fe + XFeO/DN2FeO);
    mixDiffCoeffs[gases::Fe]  = ((M_mean-XFe*M_Fe)/M_mean)/(XO2/DO2Fe + XN2/DN2Fe + XFeO/DFeFeO);
    mixDiffCoeffs[gases::FeO] = ((M_mean-XFeO*M_FeO)/M_mean)/(XO2/DO2FeO + XN2/DN2FeO + XFe/DFeFeO);
    
    return mixDiffCoeffs;
}
