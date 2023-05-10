#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

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

double M     = 21; //placeholder value for source.cpp

double A     = 6.85*1e9;           // m3/kgs
double q     = 43.28;              // must be scaled RT0/M
double Ea    = 46.37;              // must be scaled RT0
double T0    = 300.0;              // K

#endif