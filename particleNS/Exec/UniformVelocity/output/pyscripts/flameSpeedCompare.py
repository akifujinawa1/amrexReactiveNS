import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib

import os
import sys

from moviepy.editor import *
from moviepy.video.io.bindings import mplfig_to_npimage

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

colors = ['#000000', '#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D', '#f3d00d', '#9b9765']
markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'P']

# matplot subplot
fig, ax = plt.subplots(figsize=(8,8*yratio))  #fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))nrows=2,ncols=1,,dpi=100
plt.subplots_adjust(left=0.125, bottom=0.137, right=0.96, top=0.96, wspace=0.20, hspace=0.20)

condition = 2   # 1 for isobaric, 2 for isochoric
plotVar = 2     # 1 for x-t, 2 for flame speed
Nparam = 8

time = numpy.empty(((50, Nparam,2)))
location = numpy.empty(((50, Nparam,2)))
flameSpeed = numpy.empty((Nparam,2))
stDev      = numpy.empty((Nparam,2))
error        = numpy.empty(((2,Nparam,2)))
flameSpeedAvg = numpy.empty(Nparam)
concArray = [600,700,800,900,1000,1100,1200,1300]
phiArray  = numpy.empty(Nparam)

for plotTing in range(2):
    condition = plotTing+1
    for i in range(Nparam):
        concentration = (i+6)*100
        if condition == 1:
            directory = 'output/txt/1Dflame/'+str(concentration)+'/particle/'  #
        else:
            directory = 'output/txt/1DflameConfined/'+str(concentration)+'/particle/'  #
        Np = 0
        time[:,i,plotTing] = time[:,i,plotTing]*0
        location[:,i,plotTing] = location[:,i,plotTing]*0
            
        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            # checking if it is a file
            if os.path.isfile(f):
                data = np.loadtxt(f)
                # 0=time, 1=x, 2=mFe, 3=mFeO, 4=mFe3O4, 5=Tp, 6=regime
                # print(len(data[:,0]))s
                for k in range(len(data[:,0])):
                    length = len(data[:,0])
                    if (data[length-(k+1),6]-data[length-(k+2),6] == 1):
                        time[Np][i][plotTing]=(data[length-k,0])
                        location[Np][i][plotTing]=(data[length-k,1])
                        # print(time)
                        # print(location)
                        break
                Np += 1
        time[:,i,plotTing].sort()
        location[:,i,plotTing].sort()

        print(Np)
        # print(time[:,i])
        # print(location[:,i])

        interDist = (mTot0/concentration)**(1/3)
        mO2_all   = 512*dp0*interDist*interDist*Y_O2*(0.1*rhoHi+0.9*rhoLo)
        phiArray[i] = (Np*mFe0/mO2_all)/(2*M_Fe/M_O2)

        Nranges = 3
        Ndistance = 5
        offEnd = 8
        fsVals = np.empty(Nranges)
        for offset in range(Nranges):
            End = len(time[:,i,plotTing])-(offset+offEnd)
            Start = len(time[:,i,plotTing])-(offset+offEnd+Ndistance)

            # fsVals[offset] = 100*(location[End-1,i] - location[Start,i])/(time[End-1,i] - time[Start,i])

            A = np.vstack([time[Start:End,i,plotTing], np.ones(len(time[Start:End,i,plotTing]))]).T
            m, c = np.linalg.lstsq(A, location[Start:End,i,plotTing], rcond=None)[0]
            fsVals[offset] = 1e2*m
        
        flameSpeed[i,plotTing] = np.mean(fsVals)
        stDev[i,plotTing]      = np.std(fsVals)
        error[1,i,plotTing]    = np.max(fsVals)-flameSpeed[i,plotTing]
        error[0,i,plotTing]    = flameSpeed[i,plotTing]-np.min(fsVals)

        # End = len(time[:,i])-2
        # Start = len(time[:,i])-5
        # flameSpeedAvg[i] = 100*(location[End,i] - location[Start,i])/(time[End,i] - time[Start,i])
        print('Flame speed estimate for conc. = ',concentration,' is ',flameSpeed[i,plotTing],' cm/s')

        # # get individual flame speed measurements here
        # for z in range(len(location[:,1])-1):
        #     if location[z,i] == 0:
        #         flameSpeedAvg[z,i] = 0
        #     else:
        #         flameSpeedAvg[z,i] = (location[z+1,i] - location[z,i])/(time[z+1,i] - time [z,i])

        # ax.scatter(phiArray[i],flameSpeed[i],c='black',s=25,marker=markers[i],label='$'+str(concentration)+'\;\mathrm{g/cm^3}$')
        # ax.errorbar(phiArray[i],flameSpeed[i],yerr=error[:,i],fmt=markers[i],c=colors[i])

        if plotVar==1:
            ax.scatter(time[:,i,plotTing],100*location[:,i,plotTing],c=colors[i],s=15,marker=markers[i],label='$'+str(concentration)+'\;\mathrm{g/cm^3}$') #s=3
        # ax.plot(time[:,i],m*time[:,i]+c,c=colors[3],linewidth=4,label=str(concentration)) #s=3


# plt.yscale("log")
# plt.xscale("log")

# ax.ticklabel_format(useOffset=False)
if plotVar == 1:
    # ax.set_ylim([0.256,0.768])
    ax.set_ylabel(r'$x\;[\mathrm{cm}]$', fontsize=20)
    ax.set_xlabel(r'$\mathrm{time}\;[\mathrm{s}]$', fontsize=20)
if plotVar == 2:
    ax.set_ylabel(r'$\eta\;[\mathrm{cm/s}]$', fontsize=20)

    # to plot concentration on x
    ax.set_xlabel(r'$\mathrm{Concentration}\;[\mathrm{g/m^3}]$', fontsize=20)
    ax.scatter(concArray,flameSpeed[:,0],c='black',s=3,label='$\mathrm{Isobaric\;flame\;speed}$')
    ax.errorbar(concArray,flameSpeed[:,0],yerr=error[:,:,0],fmt='o',c='black',linewidth=2,label='$\sigma$')

    ax.scatter(concArray,flameSpeed[:,1],c='red',s=3,label='$\mathrm{Isochoric\;flame\;speed}$')
    ax.errorbar(concArray,flameSpeed[:,1],yerr=error[:,:,1],fmt='o',c='red',linewidth=2,label='$\sigma$')

    # to plot equivalence ratio on x
    # ax.set_xlabel(r'$\phi\;[\mathrm{-}]$', fontsize=20)
    # ax.scatter(phiArray,flameSpeed,c='black',s=3,label='$\mathrm{Reactive\;front\;speed}$')
    # ax.errorbar(phiArray,flameSpeed,yerr=error,fmt='o',c='black',linewidth=2,label='$\sigma$')

    # to plot equivalence ratio on bottom, concentration on top
    # ax.set_xlim(phiArray[0]-0.025, phiArray[8]+0.025)
    # ax.set_xlabel(r'$\phi\;[\mathrm{-}]$', fontsize=20)
    # ax.scatter(phiArray,flameSpeed,c='black',s=3,label='$\mathrm{Reactive\;front\;speed}$')

# ax2 = fig.add_axes((0.14,0.12,0.8,0.0))
# ax2.yaxis.set_visible(False) # hide the yaxis
# ax2.set_xlim(600, 1400)
# ax2.set_xticks([600,700,800,900,1000,1100,1200,1300,1400])
# ax2.set_xticklabels(['600','700','800','900','1000','1100','1200','1300','1400'])
# ax2.set_xlabel(r'$\mathrm{Concentration}\;[\mathrm{g/m^3}]$', fontsize=20)


# ax2 = ax.secondary_xaxis("top")
# ax2.set_xlabel(r'$\mathrm{Concentration}\;[\mathrm{g/m^3}]$', fontsize=20)
# concArray = [600,700,800,900,1000,1100,1200,1300,1400]

# ax2.set_xlim(600, 1400)
# ax2.set_xticks([600,700,800,900,1000,1100,1200,1300,1400])
# ax2.set_xticklabels(['600','700','800','900','1000','1100','1200','1300','1400'])



# ax.plot(concArray,flameSpeedAvg,c='red',linewidth=3,label='$\mathrm{Flame\;speed\;avg}$')

# ax.plot(location[5:end,1],flameSpeedAvg[5:end,1],c='black',linewidth=3,label='$\mathrm{Flame\;speed}$')

# print(flameSpeedAvg[:,1])


# # # ax2.set_ylim(1200,2750)
# # # ax2.set_xlim(0,0.035)

# # # # set axes limits
# # # ax2.set_xlim(0,0.01)
# # # ax[0,1].set_xlim(0,1)
# # # ax[1,0].set_xlim(0,1)
# # # ax[1,1].set_xlim(0,1)
# # # ax[0,0].set_ylim(0,1.01)
# # # ax[0,1].set_ylim(0,1.01)
# # # ax[1,0].set_ylim(0,1.01)
# # # ax[1,1].set_ylim(1.391,1.401)



# # # plt2.legend()
# # # plt2.show()

if plotVar == 1:
    ax.legend(ncol=2, loc="best", fontsize = 16)
    plt.show()
    if condition == 1:
        fig.savefig('output/plots/flame/ignition_isobaric_X-T.pdf')
    else:
        fig.savefig('output/plots/flame/ignition_isochoric_X-T.pdf')
else:
    ax.legend(ncol=1, loc="best", fontsize = 16)
    plt.show()
    if condition == 1:
        fig.savefig('output/plots/flame/flameSpeed_isobaric.pdf')
    else:
        fig.savefig('output/plots/flame/flameSpeed_isochoric.pdf')
# with open('output/txt/1Dflame/phi1/x-t.txt', 'w') as text_file:
#     for i in range(len(time)):
#         timeval = time[i]
#         locationval = location[i]
#         text_file.write(str(timeval)+' '+str(locationval)+'\n')

# text_file.close()