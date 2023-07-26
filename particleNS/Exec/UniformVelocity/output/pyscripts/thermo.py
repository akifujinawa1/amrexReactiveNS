import numpy as np
import pandas as pdf

def kO2(Tp):
    kO2_NASA = 0
    A1 = 0.77238828; A2 = 0.90875998;
    B1 = 6.9293259;  B2 = 289.86028;
    C1 = -5900.8518; C2 = -79180.433;
    D1 = 1.2202965;  D2 = 0.068622859;
    if Tp > 200 and Tp <= 1000:
        kO2_NASA = 1.0e-4*np.exp(A1*np.log(Tp) + B1/Tp + C1/Tp**2 + D1);
    else: #Tp > 1000 and Tp < 5000
        kO2_NASA = 1.0e-4*np.exp(A2*np.log(Tp) + B2/Tp + C2/Tp**2 + D2);
    return kO2_NASA

def kN2(Tp):
    kN2_NASA = 0
    A1 = 0.77238828; A2 = 0.90875998;
    B1 = 6.9293259;  B2 = 289.86028;
    C1 = -5900.8518; C2 = -79180.433;
    D1 = 1.2202965;  D2 = 0.068622859;
    if Tp > 200 and Tp <= 1000:
        kN2_NASA = 1.0e-4*np.exp(A1*np.log(Tp) + B1/Tp + C1/Tp**2 + D1);
    else: #Tp > 1000 and Tp < 5000
        kN2_NASA = 1.0e-4*np.exp(A2*np.log(Tp) + B2/Tp + C2/Tp**2 + D2);
    return kN2_NASA

def kMix(kO2, kN2, Y_O2, Y_N2):
    MO2  = 31.9988*1e-3;       # molecular mass of O2 in kg/mol
    MN2  = 28.0134*1e-3;       # molecular mass of N2 in kg/mol
    X_O2 = (Y_O2/MO2)/(Y_O2/MO2+Y_N2/MN2); # %mol fraction of oxygen in air
    X_N2  = (Y_N2/MN2)/(Y_O2/MO2+Y_N2/MN2); #  %mol fraction of other gas in air

    k_Avg = 0.5*((X_O2*kO2+X_N2*kN2) + 1/(X_O2/kO2+X_N2/kN2));
    return k_Avg

def cpO2(Tg):
    MO2 = 31.9988*1e-3;
    cp_O2=0
    M1=[31.32234,	30.03235,	20.91111]
    M2=    [-20.23531,	8.772972,	10.72071]
    M3=    [57.86644,	-3.988133,	-2.020498]
    M4=    [-36.50624,	0.788313,	0.146449]
    M5=    [-0.007374,	-0.741599,	9.245722]
    M6=    [-8.903471,	-11.32468,	5.337651]
    M7=    [246.7945,	236.1663,	237.6185]
    M8=    [0.0,	        0.0,	        0.0]

    if Tg <= 700:
        j = 1;
        cp_O2 = (1/MO2)*(M1[j] + M2[j]*(Tg/1000) + M3[j]*(Tg/1000)**2 +M4[j]*(Tg/1000)**3 + M5[j]/(Tg/1000)**2);

    elif Tg > 700 and Tg <= 2000:
        j = 2;
        cp_O2 = (1/MO2)*(M1[j] + M2[j]*(Tg/1000) + M3[j]*(Tg/1000)**2 +M4[j]*(Tg/1000)**3 + M5[j]/(Tg/1000)**2);

    else: #%Tp > 2000 and Tp <= 6000 
        j = 3;
        cp_O2 = (1/MO2)*(M1[j] + M2[j]*(Tg/1000) + M3[j]*(Tg/1000)**2 +M4[j]*(Tg/1000)**3 + M5[j]/(Tg/1000)**2);

    return cp_O2

def cpN2(Tg):
    cp_N2 = 0
    MN2 = 28.0134*1e-3;
    
    M1 = [28.98641,	19.50583,	35.51872]
    M2 =    [1.853978,	19.88705,	1.128728]
    M3 =    [-9.647459,	-8.598535,	-0.196103]
    M4 =    [16.63537,	1.369784,	0.014662]
    M5 =    [0.000117,	0.527601,	-4.553760]
    M6 =    [-8.671914,	-4.935202,	-18.97091]
    M7 =    [226.4168,	212.3900,	224.9810]
    M8 =    [0.0,	        0.0,	        0.0]

    if Tg <= 500:
        j = 1;
        cp_N2 = (1/MN2)*(M1[j] + M2[j]*(Tg/1000) + M3[j]*(Tg/1000)**2 +M4[j]*(Tg/1000)**3 + M5[j]/(Tg/1000)**2);

    elif Tg > 500 and Tg <= 2000:
        j = 2;
        cp_N2 = (1/MN2)*(M1[j] + M2[j]*(Tg/1000) + M3[j]*(Tg/1000)**2 +M4[j]*(Tg/1000)**3 + M5[j]/(Tg/1000)**2);

    else: # %Tp > 2000 && Tp <= 6000 
        j = 3;
        cp_N2 = (1/MN2)*(M1[j] + M2[j]*(Tg/1000) + M3[j]*(Tg/1000)**2 +M4[j]*(Tg/1000)**3 + M5[j]/(Tg/1000)**2);

    return cp_N2

def cpMix(cpO2,cpN2,Y_O2,Y_N2):
    cp_Mix = cpO2*Y_O2+cpN2*Y_N2
    return cp_Mix