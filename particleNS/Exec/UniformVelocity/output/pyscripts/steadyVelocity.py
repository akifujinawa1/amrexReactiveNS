import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib
from thermo import *

import os
import sys

# from moviepy.editor import *
# from moviepy.video.io.bindings import mplfig_to_npimage

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

colors = ['#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D']
color = colors

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
YO2  = M_O2*X_O2/(M_O2*X_O2 + M_N2*X_N2);   #// mass fraction of O2
YN2  = M_N2*X_N2/(M_O2*X_O2 + M_N2*X_N2);   #// mass fraction of N2
Mavg  = 1/(YO2/M_O2+YN2/M_N2);             #// average molecular weight of gas mixture
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


fig, ax = plt.subplots(figsize=(8,8*yratio))  #fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))nrows=2,ncols=1,,dpi=100
plt.subplots_adjust(left=0.16, bottom=0.16, right=0.84, top=0.9, wspace=0.20, hspace=0.20)
# plt.tight_layout()

# this solver assumes a gas density gradient via consumption of oxygen by particles,
# and that its profile is steady (1-d line)
# based on the continuity equation we can predict the 1-d velocity profile 
xlo = 0
xhi = 0.00512
Ncell = 512
dx = (xhi-xlo)/Ncell
delta = 0.0025
rhoUB = rhoLo
rhoB  = p*M_N2/(R*TpHi)
rho = np.empty(Ncell)
u = np.empty(Ncell)
for i in range(Ncell):
    x = 0.5*dx + dx*i
    if x < delta:
        a = -(rhoUB-rhoB)/(delta)
        b = rhoUB
        rho[i] = a*x+b
    else:
        rho[i] = rhoB
# rho = np.fliplr(rho) 
u[0] = 0.01
for i in range(Ncell-1):
    u[i+1] = u[i] + dx*(u[i]/rho[i])*(rho[i+1]-rho[i])/dx

xplot = np.linspace(xlo,xhi,Ncell)


ax.plot(xplot,rho,c='black',linewidth=2)
axRight = ax.twinx()
axRight.plot(xplot,u,c='red',ls='--',linewidth=2)

ax.set_ylabel(r'$\rho\;[\mathrm{kg/m^3}]$', fontsize=20)
axRight.set_ylabel(r'$u _\mathrm{g}\;[\mathrm{m/s}]$', fontsize=20)
ax.set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
# ax.errorbar(phiArray,flameSpeed,yerr=error,fmt='o',c='black',linewidth=2,label='$\sigma$')
# ax.set_xlim([phiPlotLo,phiPlotHi])
ax.legend(ncol=1, loc="best", fontsize = 14)
plt.show()