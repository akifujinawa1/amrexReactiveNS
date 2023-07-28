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

interDist = (mTot0/865)**(1/3)

nN2 = (0.1*rhoHi+0.9*rhoLo)*YN2/M_N2
 
# matplot subplot
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,8*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
plt.subplots_adjust(left=0.14, bottom=0.15, right=0.90, top=0.94, wspace=0.20, hspace=0.20)

dataAFT = np.loadtxt('output/pyscripts/aft_phi_closed.txt')

N = len(dataAFT[:,0])

phiLo = dataAFT[0,0]

# yO2offset = 0.0
# yO2lo = 0.5*((1-phiLo)*YO2+YO2)+yO2offset
# yO2hi = 0.5*(YO2)+yO2offset

# YO2first = 0.5*
# YO2next = np.empty(N-46) + 0.5*(YO2)+yO2offset


# YO2range = np.concatenate((YO2first,YO2next),axis=None)
YO2range = np.empty(N)
for iter in range(N):
    YO2range[iter] = 0.5*(YO2+max(1-dataAFT[iter,0],0)*YO2)
# print(YO2range)
# print(dataAFT[:,1])

Trange = dataAFT[:,1]
phirange = dataAFT[:,0]
eta = np.empty(N)
DToverCp     = np.empty(N)

for i in range(N):
    T = 0.5*(Trange[i]+300)    # predicted reaction zone temperature
    YO2 = YO2range[i]          # predicted reaction zone oxygen mass fraction
    Mavg  = 1/(YO2/M_O2+(1-YO2)/M_N2);    # reaction zone average molecular weight
    p = nN2*R*Trange[i]                  # adiabatic flame pressure
    D = DO2N2(T,p)

    const = (2*pi*rp0*2/(mFe0))     # constants
    DToverCp[i] = D*kMix(kO2(T),kN2(T),YO2,1-YO2)/cpMix(cpO2(T),cpN2(T),YO2,1-YO2)
    XO2term = 1/(YO2+(1-YO2)*M_O2/M_N2)

    # eta[i] = 1e2*YO2*np.sqrt(const*XO2term*DToverCp[i])   # predicted flame speed in cm/s
    eta[i] = 1e2*YO2*np.sqrt(DToverCp[i])   # predicted flame speed in cm/s

    # assuming thermal diffusion timescale >> particle combustion timescale

    # rho = p*Mavg/(R*T)
    # NpLo = 1;
    # NpHi = 20;
    # error = 1
    # tol = 1e-3
    # iteration = 0
    # Lx = 512*dp0
    # while error > tol:
    #     avg = 0.5*(NpLo+NpHi)
    #     Np = 10**(avg)
    #     mTot_all = mFeTot*Np
    #     mO2_all = mTot_all/(phirange[i]*2*M_Fe/M_O2)
    #     vol_chamber = mO2_all/(0.1*rhoHi*YO2+0.9*rhoLo*YO2)
    #     Ly = (vol_chamber/(Lx))**0.5
    #     Lz = Ly
    #     numberDensity = Np/(Lx*Ly*Lz)
    #     interDist     = (mTot0/865)**(1.0/3.0)
    #     if Ly > interDist:
    #         NpHi = avg
    #     else:
    #         NpLo = avg
    #     error = abs(Ly-interDist)/(interDist)
    #     iteration += 1
    #     if iteration > 20:
    #         print('not converging')
    #         exit()    
    # iD = 512*dp0/Np
    # volume   = 512*dp0*iD*iD
    # conc = Np*mTot0/volume
    # interDist = (mTot0/conc)**(1/3)

    # eta[i] = (kMix(kO2(T),kN2(T),YO2,1-YO2)/(rho*cpMix(cpO2(T),cpN2(T),YO2,1-YO2)))/iD


print(YO2range)


# eta = eta/max(eta)
# dataAFT = np.loadtxt('output/pyscripts/aft_phi.txt')
# ax.plot(dataAFT[:,0],dataAFT[:,1],c='red',lw=3,linestyle='--',label='$\mathrm{ad.,Fe}$-$\mathrm{to}$-$\mathrm{FeO}$')
ax.plot(phirange,eta,c='black',linewidth=3)
ax.set_ylabel(r'$\eta \;[\mathrm{cm/s}]$', fontsize=20)
ax.set_xlabel(r'$\phi \;[-]$', fontsize=20)
# ax.errorbar(phiArray,flameSpeed,yerr=error,fmt='o',c='black',linewidth=2,label='$\sigma$')
# ax.set_xlim([phiPlotLo,phiPlotHi])
ax.legend(ncol=1, loc="best", fontsize = 14)
plt.show()
# fig.savefig('output/plots/flame/postFlameTemperature.pdf')