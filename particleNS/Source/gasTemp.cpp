// #include "eulerFunc.H"
// #include "diffusionFunc.H"
// #include "AmrLevelAdv.H"
// #include "AMReX_Vector.H"
// #include "AMReX_Array.H"
// #include "AMReX_REAL.H"
// #include <AMReX_MultiFab.H>
// #include <AMReX_FArrayBox.H>
// #include <AMReX_MFIter.H>
// #include <AMReX_Amr.H>
// #include <AMReX_ParmParse.H>
// #include <AMReX_ParallelDescriptor.H>
// #include <AMReX_AmrLevel.H>

#include <cmath>
#include <iostream>
#include <algorithm> 

// Use this file to write functions required for thermodynamic and transport calculations

double hO2(const double& Tg){

    double a1,a2,a3,a4,a5,a6,a7,b1,h;
    double R = 8.3145;
        
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
    double R = 8.3145;

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


double cpO2(const double& Tg){

    double a1,a2,a3,a4,a5,a6,a7,b1,cp;
    double R = 8.3145;
        
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
    double R = 8.3145;

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
    double M_O2  = 31.9988*1e-3;       // molecular mass of O2 in kg/mol
    double M_N2  = 28.0134*1e-3;       // molecular mass of N2 in kg/mol
    double cp = Y_O2*cpO2/M_O2 + Y_N2*cpN2/M_N2;
    return cp;
}

double Hgas(const double& rho, const double& YO2, const double& YN2, const double& Tg){
    // returns gas-phase sensible enthalpy (std H of form. of O2, N2 are zero) of the 
    // binary O2-N2 gas mixture as a function of temperature, YO2, and YN2, in units
    // J/m^3 (scaled by density to make volumetric)
    double M_O2  = 31.9988*1e-3;       // molecular mass of O2 in kg/mol
    double M_N2  = 28.0134*1e-3;       // molecular mass of N2 in kg/mol
    double H = rho*(YO2*hO2(Tg)/M_O2 + YN2*hN2(Tg)/M_N2);
    return H;
}


double getMavg(const double& YO2, const double& YN2){
    double M_O2  = 31.9988*1e-3;       // molecular mass of O2 in kg/mol
    double M_N2  = 28.0134*1e-3;       // molecular mass of N2 in kg/mol
    double Mavg = 1.0/(YO2/M_O2+YN2/M_N2);
    return Mavg;
}


double Tg(const double& rho, const double& u, const double& v, \
          const double& YO2, const double& YN2, const double& ener){
    double Tg0=300,tol=1e-4,error=1;
    double dE,E,Eest,vT,Mavg,cp,Tgn,iter=0;
    double M_O2  = 31.9988*1e-3;       // molecular mass of O2 in kg/mol
    double M_N2  = 28.0134*1e-3;       // molecular mass of N2 in kg/mol

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
        Tg0 = Tgn;
        iter += 1;
        if (iter > 100){
            std::cout << "rho, u, v, YO2, YN2, ener:\n" << rho << " " << u << " " << v << " " << YO2 << \
            " " << YN2 << " " << ener << std::endl;
            std::cout << "enthalpy guess: " << enthalpy << \
            ", p guess: " << pressure << ", Tgn: " << Tgn << ", error: " << error << std::endl;
            Abort("Gas temperature iteration does not converge"); 
        }
    }

    return Tgn;

}
