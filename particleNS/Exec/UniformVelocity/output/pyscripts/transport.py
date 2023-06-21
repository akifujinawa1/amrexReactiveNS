import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib

import os
import sys

# sys.path.append(os.path.join('pyblish','plots'))
# import publish 


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

yratio = 1/1.618

# fig,ax = plt.subplots(figsize=(10,6))
fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(10,8))
fig.tight_layout(pad=2.6)


data = np.loadtxt('output/txt/transportProperties/exp.txt')

Tp = data[:,0]
kFeC = data[:,1]
kO2C = data[:,2]
kN2C = data[:,3]
kFeF = data[:,4]
kO2N = data[:,5]
kN2N = data[:,6]
muFeC = data[:,7]
muO2C = data[:,8]
muN2C = data[:,9]
muFeF = data[:,10]
muO2N = data[:,11]
muN2N  = data[:,12]
kMixC  = data[:,13]
kMixW  = data[:,14]
muMixC  = data[:,15]
muMixW  = data[:,16]
DO2C  = data[:,17]
DN2C  = data[:,18]
DFeC = data[:,19]
DFeOC  = data[:,20]
DO2CRF  = data[:,21]
DN2CRF  = data[:,22]
DFeCRF = data[:,23]
DFeOCRF  = data[:,24]

m = 100
N = int(np.floor(2800/m))
DO2red = np.empty(N)
DN2red = np.empty(N) 
DFered = np.empty(N) 
DFeOred = np.empty(N) 
TpDiff = np.empty(N) 

kO2red = np.empty(N)
kN2red = np.empty(N) 
kFered = np.empty(N) 
kMixred= np.empty(N) 
muO2red = np.empty(N)
muN2red = np.empty(N) 
muFered = np.empty(N) 
muMixred = np.empty(N) 

int = 0

for i in range(2800):
    if (np.mod(i,m) == 0):

        DO2red[int]=DO2CRF[i]
        DN2red[int]=DN2CRF[i]
        DFered[int]=DFeCRF[i]
        DFeOred[int]=DFeOCRF[i]
        TpDiff[int]=Tp[i] 
        kO2red[int]=kO2N[i]
        kN2red[int]=kN2N[i]
        kFered[int]=kFeF[i]
        kMixred[int]=kMixW[i]
        muO2red[int]=muO2N[i]
        muN2red[int]=muN2N[i]
        muFered[int]=muFeF[i]
        muMixred[int]=muMixW[i]
        int += 1 


ax[0,0].plot(Tp,kO2C,c='black',linewidth=2,linestyle='solid', label='$\mathrm{O_2,Cantera}$') 
ax[0,0].scatter(TpDiff,kO2red,c='black',s=15, label='$\mathrm{O_2,NASA}$') 

ax[0,0].plot(Tp,kN2C,c=color[0],linewidth=2,linestyle='solid', label='$\mathrm{N_2,Cantera}$') 
ax[0,0].scatter(TpDiff,kN2red,c=color[0],s=15, label='$\mathrm{N_2,NASA}$') 

ax[0,0].plot(Tp,kFeC,c='red',linewidth=2,linestyle='solid', label='$\mathrm{Fe,Cantera}$') 
ax[0,0].scatter(TpDiff,kFered,c='red',s=15, label='$\mathrm{Fe,Fitted}$') 

ax[0,1].plot(Tp,muO2C,c='black',linewidth=2,linestyle='solid', label='$\mathrm{O_2,Cantera}$') 
ax[0,1].scatter(TpDiff,muO2red,c='black',s=15, label='$\mathrm{O_2,NASA}$') 

ax[0,1].plot(Tp,muN2C,c=color[0],linewidth=2,linestyle='solid', label='$\mathrm{O_2,Cantera}$') 
ax[0,1].scatter(TpDiff,muN2red,c=color[0],s=15, label='$\mathrm{N_2,NASA}$') 

ax[0,1].plot(Tp,muFeC,c='red',linewidth=2,linestyle='solid', label='$\mathrm{Fe,Cantera}$') 
ax[0,1].scatter(TpDiff,muFered,c='red',s=15, label='$\mathrm{Fe,Fitted}$') 

ax[1,0].plot(Tp,kMixC,c='black',linewidth=2,linestyle='solid', label='$\mathrm{mix,Cantera}$') 
ax[1,0].scatter(TpDiff,kMixred,c='red',s=15, label='$\mathrm{mix,Burgoyne\;}et\;al.$') 

ax[1,1].plot(Tp,muMixC,c='black',linewidth=2,linestyle='solid', label='$\mathrm{mix,Cantera}$') 
ax[1,1].scatter(TpDiff,muMixred,c='red',s=15, label='$\mathrm{mix,Wilke}$') 

ax[0,0].set_ylabel(r'$\mathrm{k}_i\;[\mathrm{W/mK}]$', fontsize=20)
ax[0,1].set_ylabel(r'$\mu _i\;[\mathrm{kg/ms}]$', fontsize=20)

ax[1,0].set_ylabel(r'$\mathrm{k}_i\;[\mathrm{W/mK}]$', fontsize=20)
ax[1,1].set_ylabel(r'$\mu _i\;[\mathrm{kg/ms}]$', fontsize=20)

ax[1,0].set_xlabel(r'$T\;[\mathrm{K}]$', fontsize=20)
ax[1,1].set_xlabel(r'$T\;[\mathrm{K}]$', fontsize=20)

ax[0,0].ticklabel_format(axis='y',style='sci', scilimits=(-3,3))
ax[0,1].ticklabel_format(axis='y',style='sci', scilimits=(-3,3))
ax[1,0].ticklabel_format(axis='y',style='sci', scilimits=(-3,3))
ax[1,1].ticklabel_format(axis='y',style='sci', scilimits=(-3,3))

ax[0,0].set_xlim(300,3000)
ax[0,1].set_xlim(300,3000)
ax[1,0].set_xlim(300,3000)
ax[1,1].set_xlim(300,3000)
ax[0,0].set_ylim(0,0.175)
ax[0,1].set_ylim(0,1e-4)
ax[1,0].set_ylim(0,0.075)
ax[1,1].set_ylim(0,5.5e-5)
# # set axes limits
# ax2.set_xlim(0,0.01)
# ax[0,1].set_xlim(0,1)
# ax[1,0].set_xlim(0,1)
# ax[1,1].set_xlim(0,1)
# ax[0,0].set_ylim(0,1.01)
# ax[0,1].set_ylim(0,1.01)
# ax[1,0].set_ylim(0,1.01)
# ax[1,1].set_ylim(1.391,1.401)

ax[0,0].legend(loc=9, bbox_to_anchor=(0.24,1.28), fontsize="13")
# ax[0,1].legend(loc=9, bbox_to_anchor=(0.25,1.2), fontsize="12")
ax[1,0].legend(loc="best", fontsize="13")
ax[1,1].legend(loc="best", fontsize="13")

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.93, top=0.84, wspace=0.24, hspace=0.23)

plt.show()

fig.savefig('output/plots/transport/mu_k.pdf')


# plot diffusion plot here



fig2,ax2 = plt.subplots(figsize=(8,8*yratio))
# fig2.tight_layout(pad=2.6)



ax2.plot(Tp,DN2C,c='black',linewidth=2,linestyle='dashed', label='$\mathrm{N_2,Cantera}$') 
ax2.scatter(TpDiff,DN2red,c='black',linewidth=0.2, label='$\mathrm{N_2,C\;and\;E}$') 

ax2.plot(Tp,DO2C,c=color[0],linewidth=2,linestyle='solid', label='$\mathrm{O_2,Cantera}$') 
ax2.scatter(TpDiff,DO2red,c=color[0],linewidth=0.2, label='$\mathrm{O_2,C\;and\;E}$') 

ax2.plot(Tp,DFeC,c=color[1],linewidth=2,linestyle='dashdot', label='$\mathrm{Fe,Cantera}$') 
ax2.scatter(TpDiff,DFered,c=color[1],linewidth=0.2, label='$\mathrm{Fe,C\;and\;E}$') 

ax2.plot(Tp,DFeOC,c=color[2],linewidth=2,linestyle='dotted', label='$\mathrm{FeO,Cantera}$') 
ax2.scatter(TpDiff,DFeOred,c=color[2],linewidth=0.2, label='$\mathrm{FeO,C\;and\;E}$') 

# ax2.plot(Tp,DO2CRF,c='red',linewidth=2,linestyle='solid', label='$\mathrm{O_2,C\&E}$') 
# ax2.plot(Tp,DN2CRF,c='red',linewidth=2,linestyle='dashed', label='$\mathrm{N_2,C\&E}$') 
# ax2.plot(Tp,DFeCRF,c='red',linewidth=2,linestyle='dashdot', label='$\mathrm{Fe,C\&E}$') 
# ax2.plot(Tp,DFeOCRF,c='red',linewidth=2,linestyle='dotted', label='$\mathrm{FeO,C\&E}$') 



ax2.set_ylabel(r'$\mathcal{D}_i\;[\mathrm{m^2/s}]$', fontsize=20)
ax2.set_xlabel(r'$T\;[\mathrm{K}]$', fontsize=20)

ax2.ticklabel_format(axis='y',style='sci', scilimits=(-3,3))

# ax2.set_ylim(1000,3000)

# # set axes limits
# ax2.set_xlim(0,0.01)
# ax[0,1].set_xlim(0,1)
# ax[1,0].set_xlim(0,1)
# ax[1,1].set_xlim(0,1)
# ax[0,0].set_ylim(0,1.01)
# ax[0,1].set_ylim(0,1.01)
# ax[1,0].set_ylim(0,1.01)
# ax[1,1].set_ylim(1.391,1.401)

ax2.set_xlim(300,3000)
ax2.set_ylim(0,7.5e-4)
ax2.legend(loc="best", fontsize="13")

# plt2.subplots_adjust(left=0.1, bottom=0.1, right=0.93, top=0.84, wspace=0.24, hspace=0.23)

plt.show()

fig2.savefig('output/plots/transport/diff.pdf')

# plt2.legend()
# plt2.show()

