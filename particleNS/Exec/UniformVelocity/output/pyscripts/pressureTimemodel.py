import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib
from thermo import *

import os
import sys
from pathlib import Path

# from moviepy.editor import *
# from moviepy.video.io.bindings import mplfig_to_npimage

# sys.path.append(os.path.join('pyblish','plots'))
# import publish 

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
lss = ['-','--','-.']

# tfinal = 40000

# frames = 400
# fpsval = 20
# duration = frames/fpsval
# timestep = 1/duration
# framestep = tfinal/frames
# scale = framestep/timestep
 
# matplot subplot
fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(8,12*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
# fig2, ax2 = plt.subplots(nrows=3,ncols=1,figsize=(8,12*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
# fig3, ax3 = plt.subplots(nrows=2,ncols=1,figsize=(8,12*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))

plt.subplots_adjust(left=0.14, bottom=0.12, right=0.90, top=0.94, wspace=0.20, hspace=0.20)
plt.tight_layout()

Nlengths = 1
condition = 2
yIndex = 4
labelCount = 0

labels = ['$L_x=0.00512\;\mathrm{m},\;\mathrm{PRL,switch}$',
          '$L_x=0.00512\;\mathrm{m},\;k$-'r'$\beta$',
          '$L_x=0.00768\;\mathrm{m},\;\mathrm{PRL,switch}$',
          '$L_x=0.00768\;\mathrm{m},\;k$-'r'$\beta$',
          '$L_x=0.01024\;\mathrm{m},\;\mathrm{PRL,switch}$',
          '$L_x=0.01024\;\mathrm{m},\;k$-'r'$\beta$']
lsver  = ['-', '--','-.']

for Nvar in range(2):
    Nlengths = Nvar+1
    for model in range(2):
        if Nlengths == 1:
            if model == 0:
                directory =  '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft512/field/'
            else:
                directory =  '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft512kb/field/'
            TotalNp = 0
            
            for path in os.scandir(directory):
                if path.is_file():
                    TotalNp += 1
            times = np.empty(TotalNp)     # time values outputted in micros * 10
            iter = 0
            for path in os.scandir(directory):
                if path.is_file():
                    times[iter] = int(Path(path).stem)
                    iter += 1
            # time = np.linspace(0.1,59,590)
            domain = '0.00512'
            Lx = 512
            xval = 0.001275

        elif Nlengths == 2:
            if model == 0:
                directory =  '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft768/field/'
            else:
                directory =  '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft768kb/field/'
            TotalNp = 0
            for path in os.scandir(directory):
                if path.is_file():
                    TotalNp += 1
            times = np.empty(TotalNp)     # time values outputted in micros * 10
            iter = 0
            for path in os.scandir(directory):
                if path.is_file():
                    times[iter] = int(Path(path).stem)
                    iter += 1

            domain = '0.00768'
            Lx = 768
            xval = 0.001915

        elif Nlengths == 3:
            if model == 0:
                directory =  '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft1024/field/'
            else:
                directory =  '/../../../../mnt/d/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft1024kb/field/'
            TotalNp = 0
            for path in os.scandir(directory):
                if path.is_file():
                    TotalNp += 1
            times = np.empty(TotalNp)     # time values outputted in micros * 10
            iter = 0
            for path in os.scandir(directory):
                if path.is_file():
                    times[iter] = int(Path(path).stem)
                    iter += 1

            domain = '0.01024'
            Lx = 1024
            xval = 0.002555


        NT = len(times)
        times.sort()
        timePlot = np.empty(int(NT/100)+1)
        pressure = np.empty(int(NT/100)+1)
        dpdt     = np.empty(int(NT/100))
        # pressure = np.empty(NT)
        # dpdt     = np.empty(NT-1)
        i = 0
        j = 0
        while (j < NT):
        # for j in range(NT):
            timeval = str(int(times[j]))
            timePlot[i] = timeval
            data = np.loadtxt(directory+timeval+'.txt')
            # data = data[data[:, 0].argsort()]
            pressure[i] = np.mean(data[:,3])
            # print(pressure[i])
            # print(i)
            if i > 0:
                dpdt[i-1] = (pressure[i]-pressure[i-1])/(timePlot[i]-timePlot[i-1])
            i+=1
            j+=100


        # print(len(pressure))
        
        ax[0].plot(timePlot[0:len(timePlot)-1]/1e5, pressure[0:len(pressure)-1]*1e-6, color=colors[model*2],ls = lsver[Nvar], lw = 2,label=labels[labelCount])     #   label='$L_x='+domain+'\;\mathrm{m}$')
        ax[1].plot(timePlot[0:len(timePlot)-1]/1e5, dpdt*1e-6,color=colors[model*2],ls = lsver[Nvar], lw = 2) 
        labelCount+=1


ax[0].set_ylabel('$\overline{p}\;[\mathrm{MPa}]$', fontsize=20)
ax[0].legend(ncol=1, loc='lower right', fontsize = 14, frameon=False)
# ax[0].set_xlim([0,120])

ax[1].set_xlabel('$\mathrm{time}\;[\mathrm{ms}]$', fontsize=20)
ax[1].set_ylabel('$ \mathrm{d}\overline{p}/\mathrm{d}t \;[\mathrm{MPa/s}]$', fontsize=20)
# ax[1].set_xlim([0,120])


plt.show()



# if yIndex == 1:
#     ax[0].plot(data0[256:307,0]-0.00256,data0[256:307,yIndex],c='black',linewidth=2,label='$t=t_0$') 
#     ax[0].plot(data0[308:767,0]-0.00256,data0[308:767,yIndex],c='black',linewidth=2) 
#     ax[0].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2,label='$\mathrm{AFT,Fe}$-$\mathrm{to}$-$\mathrm{FeO}$') 
#     ax[0].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2,label='$\mathrm{AFT,equilibrium}$') 

#     ax[1].plot(data1[256:767,0]-0.00256,data1[256:767,yIndex],c='black',linewidth=2,label='$t=15.0\;\mathrm{ms}$') 
#     ax[1].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2) 
#     ax[1].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2)

#     ax[2].plot(data2[256:767,0]-0.00256,data2[256:767,yIndex],c='black',linewidth=2,label='$t=30.0\mathrm{ms}$') 
#     ax[2].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2) 
#     ax[2].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2) 
    
#     ax[3].plot(data3[256:767,0]-0.00256,data3[256:767,yIndex],c='black',linewidth=2,label='$t=45.0\;\mathrm{ms}$') 
#     ax[3].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2) 
#     ax[3].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2) 
    
#     ax[4].plot(data4[256:767,0]-0.00256,data4[256:767,yIndex],c='black',linewidth=2,label='$t=60.0\;\mathrm{ms}$') 
#     ax[4].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2) 
#     ax[4].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2) 
    
#     ax[0].set_ylim(0,2600)
#     ax[1].set_ylim(0,2600)
#     ax[2].set_ylim(0,2600)
#     ax[3].set_ylim(0,2600)
#     ax[4].set_ylim(0,2600)

#     ax[0].set_xlim(0,0.00512)
#     ax[1].set_xlim(0,0.00512)
#     ax[2].set_xlim(0,0.00512)
#     ax[3].set_xlim(0,0.00512)
#     ax[4].set_xlim(0,0.00512)

#     # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
#     # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
#     fig.supylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
#     ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
#     # fig.suptitle('Figure')
#     ax[0].legend(ncol=1, loc=(0.655, 0.125), fontsize = 14, frameon=False)
#     ax[1].legend(ncol=1, loc="lower left", fontsize = 14, frameon=False)
#     ax[2].legend(ncol=1, loc="lower left", fontsize = 14, frameon=False)
#     ax[3].legend(ncol=1, loc="lower left", fontsize = 14, frameon=False)
#     ax[4].legend(ncol=1, loc="lower left", fontsize = 14, frameon=False)

#     ax[0].get_xaxis().set_visible(False)
#     ax[1].get_xaxis().set_visible(False)
#     ax[2].get_xaxis().set_visible(False)
#     ax[3].get_xaxis().set_visible(False)

#     plt.show()

#     fig.savefig('output/plots/flame/thermalStructurePhi1.pdf')

# elif yIndex == 2:
#     ax[0].plot(data0[256:307,0]-0.00256,data0[256:307,yIndex],c='black',linewidth=3,label='$t=t_0$') 
#     ax[0].plot(data0[308:767,0]-0.00256,data0[308:767,yIndex],c='black',linewidth=3) 
#     ax[0].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3,label='$Y_\mathrm{O_2,0}$') 

#     ax[1].plot(data1[256:767,0]-0.00256,data1[256:767,yIndex],c='black',linewidth=3,label='$t=15.0\;\mathrm{ms}$') 
#     ax[1].plot(data1[256:767,0]-0.00256,data1[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

#     ax[2].plot(data2[256:767,0]-0.00256,data2[256:767,yIndex],c='black',linewidth=3,label='$t=30.0\mathrm{ms}$') 
#     ax[2].plot(data2[256:767,0]-0.00256,data2[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

#     ax[3].plot(data3[256:767,0]-0.00256,data3[256:767,yIndex],c='black',linewidth=3,label='$t=45.0\;\mathrm{ms}$') 
#     ax[3].plot(data3[256:767,0]-0.00256,data3[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

#     ax[4].plot(data4[256:767,0]-0.00256,data4[256:767,yIndex],c='black',linewidth=3,label='$t=60.0\;\mathrm{ms}$') 
#     ax[4].plot(data4[256:767,0]-0.00256,data4[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

#     ax[0].set_ylim(0,0.232917511457580)
#     ax[1].set_ylim(0,0.232917511457580)
#     ax[2].set_ylim(0,0.232917511457580)
#     ax[3].set_ylim(0,0.232917511457580)
#     ax[4].set_ylim(0,0.232917511457580)

#     ax[0].set_xlim(0,0.00512)
#     ax[1].set_xlim(0,0.00512)
#     ax[2].set_xlim(0,0.00512)
#     ax[3].set_xlim(0,0.00512)
#     ax[4].set_xlim(0,0.00512)

#     # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
#     # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
#     fig.supylabel(r'$Y_\mathrm{O_2}\;[\mathrm{-}]$', fontsize=20)
#     ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
#     # fig.suptitle('Figure')
#     ax[0].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
#     ax[1].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
#     ax[2].legend(ncol=1, loc="upper left", fontsize = 14, frameon=False)
#     ax[3].legend(ncol=1, loc="upper left", fontsize = 14, frameon=False)
#     ax[4].legend(ncol=1, loc="upper left", fontsize = 14, frameon=False)

#     ax[0].get_xaxis().set_visible(False)
#     ax[1].get_xaxis().set_visible(False)
#     ax[2].get_xaxis().set_visible(False)
#     ax[3].get_xaxis().set_visible(False)

#     plt.show()

#     fig.savefig('output/plots/flame/O2massFracPhi1.pdf')


# add = sum(data4[456:656,1])
# avg = add/200
# print(avg)

# print(max(data1[:,1]))