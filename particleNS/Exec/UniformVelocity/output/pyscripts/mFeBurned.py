import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib

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

print(mFe0)
 
# matplot subplot
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,8*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
plt.subplots_adjust(left=0.14, bottom=0.15, right=0.90, top=0.94, wspace=0.20, hspace=0.20)

Nparams = 9

mFeBurned = np.empty(Nparams)*0
phiArray = np.empty(Nparams)
concArray = np.empty(Nparams)

time = 30000

for i in range(Nparams):
    concentration = (i+6)*100
    concArray[i]=concentration
    Np=0
    directory = 'output/txt/1Dflame/'+str(concentration)+'/particle/'
    TotalNp = 0
    for path in os.scandir(directory):
        if path.is_file():
            TotalNp += 1

    interDist = (mTot0/concentration)**(1/3)
    mO2_all   = 512*dp0*interDist*interDist*Y_O2*(0.1*rhoHi+0.9*rhoLo)
    phiArray[i] = (TotalNp*mFe0/mO2_all)/(2*M_Fe/M_O2)

    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            data = np.loadtxt(f)
            # 0=time, 1=x, 2=mFe, 3=mFeO, 4=mFe3O4, 5=Tp, 6=regime
            # print(len(data[:,0]))s
            for k in range(len(data[:,0])):
                length = len(data[:,0])
                if (data[k,0]>time*1e-6):
                    mFe = data[k,2]
                    mFeO = data[k,3]
                    mFe3O4 = data[k,4]
                    mFeBurned[i] += mFe0-mFe
                    Np += 1
                    # print('position is:',location[Np])
                    # print('FeO mass fraction is,',YFeO[Np])
                    break

    print('conc=',concentration,', mFe burned=',mFeBurned[i])
            

ax.scatter(concArray,mFeBurned,c='black',s=15) #,label='$\mathrm{post}$-$\mathrm{flame,avg.}$')
ax.set_ylabel(r'$m_\mathrm{Fe,burned}[\mathrm{kg}]$', fontsize=20)
ax.set_xlabel(r'$\phi\;[\mathrm{-}]$', fontsize=20)
# ax.errorbar(phiArray,flameSpeed,yerr=error,fmt='o',c='black',linewidth=2,label='$\sigma$')
# ax.set_xlim([phiPlotLo,phiPlotHi])
ax.legend(ncol=1, loc="best", fontsize = 14)
plt.show()
# fig.savefig('output/plots/flame/mFeBurned.pdf')