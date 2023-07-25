import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib
import math

import os
import sys

# sys.path.append(os.path.join('pyblish','plots'))
# import publish 

yratio = 1/1.618


plt.close('all')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams["mathtext.fontset"] = 'cm'
plt.rc('font', family='serif', size='14')

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

# rcParams["axes.prop_cycle"] = cycler('color', ['#7E301E', '#B42E0F', '#FE3C0F', '#fe770f'])
colors = ['#7E301E', '#B42E0F', '#FE3C0F', '#fe770f']
color = colors

# fig,ax = plt.subplots(figsize=(10,6))






M_O2  = 31.9988*1e-3;       # molecular mass of O2 in kg/mol
M_N2  = 28.0134*1e-3;       # molecular mass of N2 in kg/mol
rhoFe     = 7874.0;         # kg/m^3                  Solid-phase iron(Fe) density
rhoFeO    = 5745.0;         # kg/m^3                  Solid-phase FeO density
rhoFe3O4  = 5170.0;         # kg/m^3                  Solid-phase Fe3O4 density
M_Fe      = 55.845*1e-3;    # kg/mol                  Molar mass of iron (Fe)
M_FeO     = 71.844*1e-3;    # kg/mol                  Molar mass of wustite (FeO)
M_Fe3O4   = 231.533*1e-3;   # kg/mol                  Molar mass of wustite (FeO)

TpHi = 1270;
TpLo = 300;
X_O2  = 0.21;               #// mole fraction of O2
X_N2  = 0.79;               #// mole fraction of N2
Y_O2  = M_O2*X_O2/(M_O2*X_O2 + M_N2*X_N2);   #// mass fraction of O2
Y_N2  = M_N2*X_N2/(M_O2*X_O2 + M_N2*X_N2);   #// mass fraction of N2
Mavg  = 1/(Y_O2/M_O2+Y_N2/M_N2);             #// average molecular weight of gas mixture
p     = 101325;
R     = 8.31446261815324;   #// J/K/mol

rhoHi = p*Mavg/(R*TpHi);
rhoLo = p*Mavg/(R*TpLo);

dp0   = 10.0e-6;  # Initial particle diameter is 10 microns
rp0   = dp0/2;
delta0 = 1.0e-3;
deltaFeO   = 0.95*delta0;
deltaFe3O4 = 0.05*delta0;
pi = 3.14159265358979

rFeO0   = rp0*(1.0-deltaFe3O4);
rFe0    = rp0*(1.0-delta0);                      # %m Initial Fe radius
mFe0    = rhoFe*(4.0/3.0)*pi*rFe0**3                 #  %kg Initial Fe mass
mFeO0   = rhoFeO*(4.0/3.0)*pi*(rFeO0**3 - rFe0**3)     #  %kg Initial FeO mass
mFe3O40 = rhoFe3O4*(4.0/3.0)*pi*(rp0**3 - rFeO0**3)     #  %kg Initial FeO mass

mTot0 = 1.0e3*(mFe0+mFeO0+mFe3O40)  # total initial particle mass in grams
mFeTot = mFe0

Lx = 0.00512   # 4 cm length of channel
Tgas = 300
rhogas = p*Mavg/(R*Tgas)

dx = dp0;

NpLo = 1;
NpHi = 20;


error = 1
tol = 1e-3
iteration = 0

conc = 1400

while error > tol:
    avg = 0.5*(NpLo+NpHi)
    Np = 10**(avg)
    print('total number of particles is: ', Np)

    # by concentration = 900-1300 g/m^3, avg 1100 g/m^3
    mTot_all = mTot0*Np    # total mass of particles in grams
    vol_chamber = mTot_all/conc
    mO2_all = vol_chamber*(0.1*rhoHi*Y_O2+0.9*rhoLo*Y_O2)

    # by equivalence ratio:

    # mTot_all = mFeTot*Np
    # mO2_all = mTot_all/(1.0*2*M_Fe/M_O2)
    # vol_chamber = mO2_all/(0.1*rhoHi*Y_O2+0.9*rhoLo*Y_O2)

    Ly = (vol_chamber/(Lx))**0.5
    Lz = Ly

    print('predicted y-domain length is: ', Ly)

    numberDensity = Np/(Lx*Ly*Lz)
    interDist     = (mTot0/865)**(1.0/3.0)

    print('predicted interparticle distance is: ', interDist)

    if Ly > interDist:
        NpHi = avg
    else:
        NpLo = avg

    error = abs(Ly-interDist)/(interDist)
    print('error: ', error)

    iteration += 1
    if iteration > 20:
        print('not converging')
        exit()

print('concentration is: ', Np*mTot0/vol_chamber)
print('Nx_cell: ', Lx/dp0)
print('Ny_cell for 2D: ', Ly/dp0)
print('Nz_cell for 3D: ', Ly/dp0)

print('number of particles: ', Np)
Np = math.ceil(Np)

spacing_x = math.floor((Lx/dp0)/Np)

print('cells between particles in x: ', spacing_x)

# with open('setupScripts/locations.txt', 'w') as text_file:
#     for i in range(math.floor(Np)):
#         x = 256+math.ceil(0.5*spacing_x)+spacing_x*i
#         text_file.write(repr(x)+'\n')

# text_file.close()

equiv = (Np*mFe0/mO2_all)/(2*M_Fe/M_O2)
print('equivalence ratio of sim: ', equiv)

# yratio = 1/1.618
# data = np.loadtxt('setupScripts/locations.txt')
# fig,ax = plt.subplots(figsize=(10,2))

# ax.scatter(data[:,0],data[:,1],c='black',linewidth=1) 
# plt.show()