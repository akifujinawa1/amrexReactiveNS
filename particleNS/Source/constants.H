#ifndef CONSTANTS_H
#define CONSTANTS_H

// #include <cmath>
#include <math.h>

// Gas-phase constants

double Gamma = 1.4;
double R     = 8.31446261815324;   // J/K/mol
double Pr    = 0.75;               // dimensionless
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
double Dpre  = (1e-7*pow(1/M_O2 + 1/M_N2,0.5))/pow(pow(vO2,(1.0/3.0))+pow(vN2,(1.0/3.0)),2);


// Particle constants

double delta0    = 1e-3;           // Initial oxide layer thickness to particle size ratio
double rhoFe     = 7874.0;         // kg/m^3                  Solid-phase iron(Fe) density
double rhoFeO    = 5745.0;         // kg/m^3                  Solid-phase FeO density
double rhoFe3O4  = 5170.0;         // kg/m^3                  Solid-phase Fe3O4 density
double M_Fe      = 55.845*1e-3;    // kg/mol                  Molar mass of iron (Fe)
double M_FeO     = 71.844*1e-3;    // kg/mol                  Molar mass of wustite (FeO)
double M_Fe3O4   = 231.533*1e-3;   // kg/mol                  Molar mass of wustite (FeO)
double HfuFe     = 209*1e3;           // J/kg                 Heat of fusion to melt Fe
double HfuFeO    = 7300*4.184/M_FeO;  // J/kg                 Heat of fusion to melt FeO
double HfuFe3O4  = 596*1e3;           // J/kg                 Heat of fusion to melt Fe3O4

double k0FeOs    = 2.670e-4;       // m^2/s                   Pre-exponential factor for kinetic rate of FeO oxide layer growth
double k0Fe3O4s  = 1.027e-6;       // m^2/s                   Pre-exponential factor for kinetic rate of Fe3O4 oxide layer growth
double TaFeOs    = 20319;          // K                       Activation temperature for kinetic rate of FeO oxide layer growth
double TaFe3O4s  = 21310;          // K                       Activation temperature for kinetic rate of Fe3O4 oxide layer growth

double nFeFeO = M_Fe/M_FeO;
double nFeFe3O4 = 3.0*M_Fe/M_Fe3O4;
double nO2FeO = 0.5*M_O2/M_FeO;
double nO2Fe3O4 = 2.0*M_O2/M_Fe3O4;   

double pi    = M_PI;

double M     = 21; //placeholder value for source.cpp

double A     = 6.85*1e9;           // m3/kgs
double q     = 43.28;              // must be scaled RT0/M
double Ea    = 46.37;              // must be scaled RT0
double T0    = 300.0;              // K

#endif