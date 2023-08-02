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

# mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

colors = ['#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D']
color = colors
markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'P']


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
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.90, top=0.94, wspace=0.20, hspace=0.20)

condition = 1
if condition == 1:
    concentrations = [725,750,800,900,1000,1100,1200,1300,1400,1500,1600]
elif condition == 2:
    concentrations = [600,700,800,900,1000,1100,1200,1300,1400]

Nparams = len(concentrations)

mFeBurned = np.empty(Nparams)
totalCp = np.empty(Nparams)
phiArray = np.empty(Nparams)
concArray = np.empty(Nparams)

if condition == 1:
    folder = '1Dflame'
else:
    folder = '1DflameConfined'

time = 20000

for i in range(Nparams):
    concentration = concentrations[i]
    concArray[i]=concentration
    Np=0
    totalmFe = 0
    directory = 'output/txt/'+folder+'/'+str(concentration)+'/particle/'
    TotalNp = 0
    for path in os.scandir(directory):
        if path.is_file():
            TotalNp += 1

    fixedID    = (mTot0/865)**(1/3)
    cpParticle = 0.1*TotalNp*1.856370123362867e-09 + 0.9*TotalNp*3.091198097151238e-09
    cpGas      = (0.1*rhoHi*cpMix(cpO2(1270),cpN2(1270),Y_O2,Y_N2)+0.9*rhoLo*cpMix(cpO2(300),cpN2(300),Y_O2,Y_N2))
    fixedVol   = fixedID*fixedID*dp0*512
    totalCp[i] = cpParticle+cpGas*fixedVol

    interDist = (mTot0/865)**(1/3)
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
                length = len(data[:,0])-1
                if (data[length,0] < time*1e-6):
                    mFe = data[length,2]
                    mFeO = data[length,3]
                    mFe3O4 = data[length,4]
                    totalmFe = totalmFe + mFe0-mFe
                    Np += 1
                    break
                if (data[k,0]>time*1e-6):
                    mFe = data[k,2]
                    mFeO = data[k,3]
                    mFe3O4 = data[k,4]
                    totalmFe = totalmFe + mFe0-mFe
                    
                    Np += 1
                    # print('position is:',location[Np])
                    # print('FeO mass fraction is,',YFeO[Np])
                    break
                
    mFeBurned[i] = totalmFe
    print('conc=',concentration,', mFe burned=',mFeBurned[i])

# totalmassO2 = phiArray*0+mO2_all*2*M_Fe/M_O2
# print(mO2_all*2*M_Fe/M_O2)
# ax.plot(phiArray,totalmassO2,c='black',linewidth=3,linestyle='dashed',label='$m_\mathrm{Fe,st.}$')
ax.scatter(phiArray,mFeBurned,c='black',marker=markers[0],s=20,label='$m_\mathrm{Fe}$')
# ax.set_ylim([4.2e-11,9.8e-11])
ax.set_ylabel(r'$m_\mathrm{Fe,burned}\;[\mathrm{kg}]$', fontsize=20)
ax.set_xlabel(r'$\phi\;[\mathrm{-}]$', fontsize=20)

ax2 = ax.twinx()
ax2.scatter(phiArray,totalCp,c='red',marker=markers[1],s=20,label='$c_{p,\mathrm{total}}$')
ax2.set_ylabel(r'$c_{p,\mathrm{total}}\;[\mathrm{J/K}]$', fontsize=20)
# ax.set_xlabel(r'$\phi\;[\mathrm{-}]$', fontsize=20)

# ax.errorbar(phiArray,flameSpeed,yerr=error,fmt='o',c='black',linewidth=2,label='$\sigma$')
# ax.set_xlim([phiPlotLo,phiPlotHi])
ax.legend(ncol=1, loc="best", fontsize = 16, frameon=False)
ax2.legend(ncol=1, loc=(0.016,0.78), fontsize = 16, frameon=False)
plt.show()

if condition == 1:
    fig.savefig('output/plots/flame/burnMassVsCp_isobaric.pdf')
else:
    fig.savefig('output/plots/flame/burnMassVsCp_isochoric.pdf')