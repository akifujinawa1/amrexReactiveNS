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
ignDelayStart = 1e-14
dTdtStart = 5e4


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


        
totalN = 0;
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        totalN+=1

time = numpy.empty(totalN)
location = numpy.empty(totalN)
startTime = np.empty(totalN)
ignDelay = np.empty(totalN)

timekb = numpy.empty(totalN)
locationkb = numpy.empty(totalN)
startTimekb = np.empty(totalN)
ignDelaykb = np.empty(totalN)

mFeval   = numpy.empty(totalN) 
dTdt     = numpy.empty(totalN)
Tign     = numpy.empty(totalN)
files = ['']*totalN
    
Np = 0

        
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    print(f)
    files[Np] = filename
    # print(files[Np])
    # print(f)
    
    # checking if it is a file
    if os.path.isfile(f):
        data = np.loadtxt(f)
        # 0=time, 1=x, 2=mFe, 3=mFeO, 4=mFe3O4, 5=Tp, 6=regime
        # print(len(data[:,0]))s
        startCheck = 0
        endCheck   = 0
        
        for k in range(len(data[:,0])-2):
            dmFeOdt = (data[k+1,3]-data[k,3])/(data[k+1,0]-data[k,0])
            dTdt0   = (data[k+1,5]-data[k,5])/(data[k+1,0]-data[k,0])
            
            # if (dmFeOdt > ignDelayStart) and (startCheck == 0):
            #     # print(dmFeOdt)
            #     startTime[Np] = data[k,0]
            #     startCheck = 1
            if (dTdt0 > dTdtStart) and (startCheck == 0):
                startTime[Np] = data[k,0]
                startCheck = 1
            length = len(data[:,0])-1
            if (data[length-(k),6]-data[length-(k+1),6] == 1) and (endCheck == 0):
                time[Np]=(data[length-k,0])
                location[Np]=(data[length-k,1])
                mFeval[Np]=(data[length-k,2])
                dTdt[Np]  = (data[length-k,5]-data[(length-k)-1,5])/(data[length-k,0]-data[(length-k)-1,0])
                # print(dTdt[Np])
                Tign[Np] = data[length-k,5]
                # print(location)
                endCheck = 1
            # if (data[k,5] > 1870) and (endCheck == 0):
            #     time[Np]=(data[k,0])
            #     location[Np]=(data[k,1])
            #     mFeval[Np]=(data[k,2])
            #     # print(time)
            #     # print(location)
            #     endCheck = 1
            if (startCheck == 1) and (endCheck == 1):
                break
        Np += 1

if domainL == 0:
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft512kb/particle/'
elif domainL == 1:
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft768kb/particle/'
elif domainL == 2:
    directory = '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft1024kb/particle/'

# Np = 0

# print(startTime)

for i in range(Np):
    # print(files[Np].decode(encoding))
    f = os.path.join(directory, files[i])
    # f =  files[i]
    print(f)
    data = np.loadtxt(f)
    # print(mFeval[i])
    startCheck = 0
    endCheck = 0
    # print(dTdt[i])
    dTdtkb   = numpy.empty(len(data[:,0]))
        
    for k in range(len(data[:,0])-2):
        dmFeOdt = (data[k+1,3]-data[k,3])/(data[k+1,0]-data[k,0])
        dTdt0   = (data[k+1,5]-data[k,5])/(data[k+1,0]-data[k,0])
        # if (dmFeOdt > ignDelayStart) and (startCheck == 0):
        #     startTimekb[i] = data[k,0]
        #     startCheck = 1
        dTdtkb[k] = dTdt0
        if (dTdt0 > dTdtStart) and (startCheck == 0):
            startTimekb[i] = data[k,0]
            startCheck = 1  
            # break
        # if (data[k,2]<mFeval[i]) and (endCheck == 0):
        #     timekb[i]=(data[k,0])
        #     locationkb[i]=(data[k,1])
        #     endCheck = 1
        # if (data[k,5] > 1870) and (endCheck == 0):
        #     timekb[i]=(data[k,0])
        #     locationkb[i]=(data[k,1])
        #     endCheck = 1
        # dTdt1 = (data[k+1,5]-data[k,5])/(data[k+1,0]-data[k,0])
        # print(dTdt0)
        # if (dTdt0 > dTdt[i]) and (endCheck == 0):
        #     # print(dTdt0)
        # if (dTdtkb[k] < dTdtkb[k-1]) and (endCheck == 0):
        #     timekb[i]=(data[k,0])
        #     locationkb[i]=(data[k,1])
        #     endCheck = 1
        if (data[k,5] > Tign[i]) and (endCheck == 0):
            timekb[i]=(data[k,0])
            locationkb[i]=(data[k,1])
            endCheck = 1
        if (startCheck == 1) and (endCheck == 1):
            break
    # Np += 1
    # print(np.max(dTdtkb[0:k]))
    # kval = np.array(dTdtkb).argmax()
    # timekb[i]=(data[kval,0])
    # locationkb[i]=(data[kval,1])

    
# time1 = time[0:Np,i]
# location1 = location[0:Np,i]
# time1.sort()
# location1.sort()
time.sort()
location.sort()
startTime.sort()

timekb.sort()
locationkb.sort()
startTimekb.sort()

for i in range(Np):
    ignDelay[i] = time[i] - startTime[i]
    ignDelaykb[i] = timekb[i] - startTimekb[i]
    print('PRL:',ignDelay[i],'k-beta:',ignDelaykb[i])

print('Average ignition delay time for PRL is ',np.mean(ignDelay),'s')
print('Average ignition delay time for k-beta is ',np.mean(ignDelaykb),'s')

# ax.scatter(time,ignDelay,c=colors[0],s=15,marker=markers[0],label=labels[0]) #s=3
# ax.scatter(timekb,ignDelay,c=colors[2],s=15,marker=markers[2],label=labels[1]) #s=3
# ax.set_xlabel(r'$\mathrm{time}\;[\mathrm{s}]$', fontsize=20)

# ax.plot(location,time,c=colors[0],lw=2,marker=markers[0],label=labels[0]) #s=3
# ax.plot(locationkb,timekb,c=colors[2],lw=2,marker=markers[2],label=labels[1]) #s=3

ax.scatter(location,ignDelay,c=colors[0],s=15,marker=markers[0],label=r'$\tau _\mathrm{ign.},\;\mathrm{PRL\;(switch)}$') #s=3
ax.scatter(locationkb,ignDelaykb,c=colors[2],s=15,marker=markers[2],label=r'$\tau _\mathrm{ign.},\;k$-$\beta$') #s=3
ax.set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)


if domainL == 0:
    # ax.set_ylim([0,1.1*np.max(ignDelay)])
    ax.set_xlim([0,0.00512])
# elif domainL == 1:
#     # ax.set_ylim([0,0.768])
#     # ax.set_xlim([0,0.07])
# elif domainL == 2:
#     # ax.set_ylim([0,1.024])
#     # ax.set_xlim([0,0.09])

ax.set_ylabel(r'$\tau _\mathrm{ign.}\;[\mathrm{s}]$', fontsize=20)
ax.legend(ncol=1, loc="best", fontsize = 16, frameon = False )
plt.show()