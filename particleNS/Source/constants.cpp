
#include "constants.H"
#include <math.h>
#include <cmath>

// Gas-phase constants

double Gamma = 1.4;
double R     = 8.31446261815324;   // J/K/mol
double one_atm_Pa = 101325.0;      // pa

double M_O2  = 31.9988*1e-3;       // molecular mass of O2 in kg/mol
double M_N2  = 28.0134*1e-3;       // molecular mass of N2 in kg/mol
double X_O2  = 0.21;               // mole fraction of O2
double X_N2  = 0.79;               // mole fraction of N2
double Y_O2  = M_O2*X_O2/(M_O2*X_O2 + M_N2*X_N2);   // mass fraction of O2
double Y_N2  = M_N2*X_N2/(M_O2*X_O2 + M_N2*X_N2);   // mass fraction of N2
double Mavg  = 1/(Y_O2/M_O2+Y_N2/M_N2);             // average molecular weight of gas mixture

double vO2   = 16.6;
double vN2   = 17.9;
double Dpre  = (1e-7*pow(1.0/M_O2 + 1.0/M_N2,0.5))/pow(pow(vO2,(1.0/3.0))+pow(vN2,(1.0/3.0)),2.0);


// Particle constants

double delta0    = 1e-3;           // Initial oxide layer thickness to particle size ratio
double rhoFe     = 7874.0;         // kg/m^3                  Solid-phase iron(Fe) density
double rhoFeO    = 5745.0;         // kg/m^3                  Solid-phase FeO density
double rhoFe3O4  = 5170.0;         // kg/m^3                  Solid-phase Fe3O4 density
double M_Fe      = 55.845*1e-3;    // kg/mol                  Molar mass of iron (Fe)
double M_FeO     = 71.844*1e-3;    // kg/mol                  Molar mass of wustite (FeO)
double M_Fe3O4   = 231.533*1e-3;   // kg/mol                  Molar mass of magnetite (Fe3O4)
double M_Fe2O3   = 159.688*1e-3;   // kg/mol                  Molar mass of hematite (Fe2O3)
double HfuFe     = 209*1e3;           // J/kg                 Heat of fusion to melt Fe
double HfuFeO    = 7300*4.184/M_FeO;  // J/kg                 Heat of fusion to melt FeO
double HfuFe3O4  = 596*1e3;           // J/kg                 Heat of fusion to melt Fe3O4

double k0FeOs    = 2.670e-4;       // m^2/s                   Pre-exponential factor for kinetic rate of FeO oxide layer growth
double k0Fe3O4s  = 1.027e-6;       // m^2/s                   Pre-exponential factor for kinetic rate of Fe3O4 oxide layer growth
double TaFeOs    = 20319;          // K                       Activation temperature for kinetic rate of FeO oxide layer growth
double TaFe3O4s  = 21310;          // K                       Activation temperature for kinetic rate of Fe3O4 oxide layer growth

double qFeOs     = 3.787e6;        // J/kg
double qFeOl     = 3.473e6;        // J/kg
double qFe3O4s   = 4.841e6;        // J/kg
double HformFe2O3s = -825.50;      // J/mol
double qFe2O3s   = -HformFe2O3s/M_Fe2O3; //J/kg
double HformFeOg = -251.04;
double qFeOg     = -HformFeOg/M_FeO; //J/kg

double pi    = M_PI;


double nFeFeO   = M_Fe/M_FeO;
double nFeFe3O4 = 3.0*M_Fe/M_Fe3O4;
double nO2FeO   = 0.5*M_O2/M_FeO;
double nO2Fe3O4 = 2.0*M_O2/M_Fe3O4;   

double rp0        = 0.5*dp0;
double deltaFeO   = 0.95*delta0;
double deltaFe3O4 = 0.05*delta0;

double rFeO0      = rp0*(1-deltaFe3O4);
double rFe0       = rp0*(1-delta0);
double mFe0       = rhoFe*(4.0/3.0)*pi*pow(rFe0,3.0);
double mFeO0      = rhoFeO*(4.0/3.0)*pi*(pow(rFeO0,3.0)-pow(rFe0,3.0));
double mFe3O40    = rhoFe3O4*(4.0/3.0)*pi*(pow(rp0,3.0)-pow(rFeO0,3.0));
double interDist  = pow(1.0e3*(mFe0+mFeO0+mFe3O40)/1100.0,1.0/3.0);


// Constants for multi-species diffusion of O2, N2, Fe, FeO
double kb      = 1.380649e-23;
double sigmaO2 = 3.46;
double sigmaN2 = 3.621; 
double sigmaFe = 4.3;
double sigmaFeO = 4.3;
double epsO2    = 107.40*kb;
double epsN2    = 97.53*kb;
double epsFe    = 3000.0*kb;
double epsFeO   = 3000.0*kb;

double sigmaO2N2  = 0.5*(sigmaO2+sigmaN2);
double sigmaO2Fe  = 0.5*(sigmaO2+sigmaFe);
double sigmaO2FeO = 0.5*(sigmaO2+sigmaFeO);
double sigmaN2Fe  = 0.5*(sigmaFe+sigmaN2);
double sigmaN2FeO = 0.5*(sigmaFeO+sigmaN2);
double sigmaFeFeO = 0.5*(sigmaFe+sigmaFeO);

double epsO2N2  = sqrt(epsO2*epsN2);
double epsO2Fe  = sqrt(epsO2*epsFe);
double epsO2FeO = sqrt(epsO2*epsFeO);
double epsN2Fe  = sqrt(epsFe*epsN2);
double epsN2FeO = sqrt(epsFeO*epsN2);
double epsFeFeO = sqrt(epsFe*epsFeO);

double M_O2N2   = 2.0/(1.0/(1e3*M_O2) + 1.0/(1e3*M_N2));   // need these quantities in g/mol for diffusion calculations
double M_O2Fe   = 2.0/(1.0/(1e3*M_O2) + 1.0/(1e3*M_Fe));
double M_O2FeO  = 2.0/(1.0/(1e3*M_O2) + 1.0/(1e3*M_FeO));
double M_N2Fe   = 2.0/(1.0/(1e3*M_Fe) + 1.0/(1e3*M_N2));
double M_FeFeO  = 2.0/(1.0/(1e3*M_Fe) + 1.0/(1e3*M_FeO));
double M_N2FeO  = 2.0/(1.0/(1e3*M_N2) + 1.0/(1e3*M_FeO));


double M     = 21; //placeholder value for source.cpp
double A     = 6.85*1e9;           // m3/kgs
double q     = 43.28;              // must be scaled RT0/M
double Ea    = 46.37;              // must be scaled RT0
double T0    = 300.0;              // K
