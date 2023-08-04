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

# mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

colors = ['#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D']
color = colors

# tfinal = 40000

# frames = 400
# fpsval = 20
# duration = frames/fpsval
# timestep = 1/duration
# framestep = tfinal/frames
# scale = framestep/timestep
 
# matplot subplot
fig1, ax1 = plt.subplots(nrows=4,ncols=1,figsize=(8,16*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
plt.subplots_adjust(left=0.14, bottom=0.08, right=0.90, top=0.97, wspace=0.20, hspace=0.20)

condition = 1
setConcentrations = 1

for i in range(4):
    if setConcentrations == 1:
        concentrations = [600,700,800,900]
        concentration = concentrations[i]
    else:
        concentrations = [1000,1100,1200,1300]
        concentration = concentrations[i]
    time = 59000
    #27600

    if condition == 1:
        folder = '1Dflame'
    else:
        folder = '1DflameConfined'

    data4 = np.loadtxt('output/txt/'+folder+'/'+str(concentration)+'/field/'+str(time)+'.txt')
    data4 = data4[data4[:, 0].argsort()]

    directory = 'output/txt/'+folder+'/'+str(concentration)+'/particle/'

    TotalNp = 0

    for path in os.scandir(directory):
        if path.is_file():
            TotalNp += 1

    location = np.empty(TotalNp)
    YFe      = np.empty(TotalNp)
    YFeO     = np.empty(TotalNp)

    matrix = np.empty((TotalNp,2))

    Np = 0

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
                    # location[Np]=data[k,1]
                    matrix[Np,0] = data[k,1] 
                    mFe = data[k,2]
                    mFeO = data[k,3]
                    mFe3O4 = data[k,4]
                    # YFe[Np] = mFe/(mFe+mFeO+mFe3O4)
                    # YFeO[Np] = mFeO/(mFe+mFeO+mFe3O4)
                    matrix[Np,1] = mFeO/(mFe+mFeO+mFe3O4)
                    Np += 1
                    # print('position is:',location[Np])
                    # print('FeO mass fraction is,',YFeO[Np])
                    break
            

    matrix = matrix[matrix[:,0].argsort()]
    # YFeO = [x for _,x in sorted(zip(YFeO,location))]
    location = matrix[:,0]
    YFeO     = matrix[:,1]
    # for olp in range(len(matrix[:,1])-1):
    #     if (YFeO[olp] < 1e-4):
    #         YFeO[olp] == 0
    #     elif YFeO[olp] > 2:
    #         YFeO[olp] == 1
    # location.sort()
    # sorted(YFeO, reverse = True)
    if condition == 1 and i == 0:
        YFeO[15] = 0

    if condition == 1:
        start = 256
        end = 767
        offset = 0.256
    else:
        start = 0
        end = 511
        offset = 0

    location = location*1e2 - offset
    # print(location)
    print(YFeO)
    # print(Np)
    # print(data4[:,2])

    interDist = (mTot0/865)**(1/3)
    mO2_all   = 512*dp0*interDist*interDist*Y_O2*(0.1*rhoHi+0.9*rhoLo)
    phival    = (TotalNp*mFe0/mO2_all)/(2*M_Fe/M_O2)
    phival    = round(phival,2)



    if i == 0:
        ax1[i].plot(data4[start:end,0]*1e2-offset,data4[start:end,2],c='black',linewidth=3,label='$Y_\mathrm{O_2}$') 
        if setConcentrations == 1:
            ax1[i].legend(ncol=1, loc=(0,0.15), fontsize = 16, frameon=False)
        else: 
            ax1[i].legend(ncol=1, loc=(0,0.23), fontsize = 16, frameon=False)    
    else:
        ax1[i].plot(data4[start:end,0]*1e2-offset,data4[start:end,2],c='black',linewidth=3)
        

    # ax[i].set_ylabel(r'$Y_\mathrm{O_2}\;[-]$', fontsize=20)
    ax1[i].set_xlabel(r'$x\;[\mathrm{cm}]$', fontsize=20)
    ax1[i].set_ylim(-0.00235,0.235)
    ax1[i].set_xlim([0,0.512])
    

    textstr = '$'+str(concentration)+'\;\mathrm{g/m^3},\;\phi \sim'+str(phival)+'$'

    if setConcentrations == 1:
        if i == 3:
            ax1[i].text(0.005, 0.04, textstr, fontsize=16)
        else:
            ax1[i].text(0.005, 0.01, textstr, fontsize=16)
    else:
        if i == 0:
            ax1[i].text(0.005, 0.03, textstr, fontsize=16)
        elif i == 1:
            ax1[i].text(0.005, 0.03, textstr, fontsize=16)
        else:
            ax1[i].text(0.005, 0.03, textstr, fontsize=16)

    # ax2.scatter(location, YFeO, color = colors[1], s=15) #s=15  linewidth=3
    axR = ax1[i].twinx()
    if i == 0:
        axR.plot(location, YFeO, '--o', color=colors[2], label='$Y_\mathrm{FeO}$') #s=15  linewidth=3
        if setConcentrations == 1:
            axR.legend(ncol=1, loc=(0.2,0.17), fontsize = 16, frameon=False)
        else:
            axR.legend(ncol=1, loc=(0.2,0.24), fontsize = 16, frameon=False)
    else:
        if (i == 3):
            axR.plot(location, YFeO, '--o', color=colors[2], linewidth = 1)
            # axR.scatter(location[len(location)-1], YFeO[len(location)-1],marker= 'o', color=colors[2])
        else: 
            axR.plot(location, YFeO, '--o', color=colors[2], linewidth = 1)
    axR.set_ylim([-0.01,1.03])
    # axR.set_xlim([0,0.00512])
    # ax2.tick_params(axis ='y', labelcolor = colors[1])

    if i == 1:
        ax1[i].set_ylabel(r'$Y_\mathrm{O_2}\;[-]$',y = -0.04, fontsize=20)
        axR.set_ylabel(r'$Y_\mathrm{FeO}\;[-]$',y = -0.04, color=colors[2], fontsize=20)
    
# fig.supylabel(r'$Y_\mathrm{O_2}\;[-]$', fontsize=20)
# fig2 = fig.twinx()
# fig2.supylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
ax1[0].get_xaxis().set_visible(False)
ax1[1].get_xaxis().set_visible(False)
ax1[2].get_xaxis().set_visible(False)


plt.show()

if setConcentrations == 1:
    if condition == 1:
        fig1.savefig('output/plots/flame/yO2yFeO_lean_isobaric.pdf')
    else:
        fig1.savefig('output/plots/flame/yO2yFeO_lean_isochoric.pdf')
else:
    if condition == 1:
        fig1.savefig('output/plots/flame/yO2yFeO_rich_isobaric.pdf')
    else:
        fig1.savefig('output/plots/flame/yO2yFeO_rich_isochoric.pdf')