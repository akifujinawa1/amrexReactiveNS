import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib
from thermo import *

import os
import sys

from moviepy.editor import *
from moviepy.video.io.bindings import mplfig_to_npimage

# sys.path.append(os.path.join('pyblish','plots'))
# import publish 

# CHOOSE PLOTTING PARAMETERS HERE

condition = 2   # 1 for isobaric, 2 for isochoric
domainL = 0
particle = '25-0.txt'


yratio = 1/1.618

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


plt.close('all')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams["mathtext.fontset"] = 'cm'
plt.rc('font', family='serif', size='14')

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False


colors = ['#000000', '#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D', '#f3d00d', '#9b9765', '#000000', '#9B0909', '#000000']
markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'P', 'o', 'v', 'o']

# matplot subplot


fig, ax = plt.subplots(figsize=(8,8*yratio))  #fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))nrows=2,ncols=1,,dpi=100
plt.subplots_adjust(left=0.17, bottom=0.137, right=0.91, top=0.9, wspace=0.20, hspace=0.20)




if domainL == 0:
    labels = ['$L_x=0.00512\;\mathrm{m},\;\mathrm{PRL\;(switch)}$',
              '$L_x=0.00512\;\mathrm{m},\;k$-'r'$\beta$']
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft512/particle/'
elif domainL == 1:
    labels = ['$L_x=0.00768\;\mathrm{m},\;\mathrm{PRL\;(switch)}$',
              '$L_x=0.00768\;\mathrm{m},\;k$-'r'$\beta$']
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft768/particle/'
elif domainL == 2:
    labels = ['$L_x=0.01024\;\mathrm{m},\;\mathrm{PRL\;(switch)}$',
              '$L_x=0.01024\;\mathrm{m},\;k$-'r'$\beta$']
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft1024/particle/'



data = np.loadtxt(directory+particle)

# 0=time, 1=x, 2=mFe, 3=mFeO, 4=mFe3O4, 5=Tp, 6=regime

ax.plot(data[:,0],data[:,5],c=colors[0],lw=2,label=r'$\mathrm{PRL\;(switch)}$') #s=3

if domainL == 0:
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft512kb/particle/'
elif domainL == 1:
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft768kb/particle/'
elif domainL == 2:
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft1024kb/particle/'

data = np.loadtxt(directory+particle)

# 0=time, 1=x, 2=mFe, 3=mFeO, 4=mFe3O4, 5=Tp, 6=regime

ax.plot(data[:,0],data[:,5],c=colors[2],lw=2,label=r'$k$-$\beta$') #s=3

ax.set_xlabel(r'$\mathrm{time}\;[\mathrm{s}]$', fontsize=20)


# if domainL == 0:
#     ax.set_ylim([0,1.1*np.max(ignDelay)])
#     ax.set_xlim([0,0.00512])
# elif domainL == 1:
#     # ax.set_ylim([0,0.768])
#     # ax.set_xlim([0,0.07])
# elif domainL == 2:
#     # ax.set_ylim([0,1.024])
#     # ax.set_xlim([0,0.09])

ax.set_ylabel(r'$T_\mathrm{p}\;[\mathrm{K}]$', fontsize=20)
ax.legend(ncol=1, loc="best", fontsize = 16, frameon = False )
plt.show()