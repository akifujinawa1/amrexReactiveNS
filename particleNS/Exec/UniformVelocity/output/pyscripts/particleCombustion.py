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

yratio = 1/1.618


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

# fig,ax = plt.subplots(figsize=(10,6))
fig2,ax2 = plt.subplots(figsize=(8,8*yratio))



dataRef  = np.loadtxt('output/txt/particleCombustion/dataZeroD.txt')
data  = np.loadtxt('output/txt/particleCombustion/data1270.txt')
dataNing  = np.loadtxt('output/txt/particleCombustion/dataNingTp.txt')
dataNingLo  = np.loadtxt('output/txt/particleCombustion/dataNingTpLo.txt')
dataNingHi  = np.loadtxt('output/txt/particleCombustion/dataNingTpHi.txt')


ax2.plot(data[:,0],data[:,1],c='red',linewidth=1,label='$T_\mathrm{p,0}=1270\;\mathrm{K}$') 
ax2.plot(dataRef[:,0],dataRef[:,1],c='black',linewidth=1,label='$\mathrm{Ref.}$') 
ax2.plot(1e-3*(dataNing[2:80,0]+2.4834281),dataNing[2:80,1],c='black',linestyle='dotted',linewidth=1,label='$\mathrm{Exp.,\;Ning\;}et\;al.}$') 
# ax2.fill_between(1e-3*(dataNingLo[0:80,0]+3.218259281),dataNingLo[0:80,1],dataNingHi[0:80,1]) 

# # set labels
# ax[1,0].set_xlabel(r'$x$')
# ax[1,1].set_xlabel(r'$x$')
# ax[0,0].set_ylabel(r'$\rho$')
# ax[0,1].set_ylabel(r'$u$')
# ax[1,0].set_ylabel(r'$p$')
# ax[1,1].set_ylabel(r'$\gamma$')
# ax.set_ylabel(r'$Y_\mathrm{O_2}$', fontsize=20)
# ax.set_xlabel(r'$x$', fontsize=20)
ax2.set_ylabel(r'$T_\mathrm{p}\;[\mathrm{K}]$', fontsize=20)
ax2.set_xlabel(r'$\mathrm{time}\;[\mathrm{s}]$', fontsize=20)

ax2.set_ylim(1200,2750)
ax2.set_xlim(0,0.035)

# # set axes limits
# ax2.set_xlim(0,0.01)
# ax[0,1].set_xlim(0,1)
# ax[1,0].set_xlim(0,1)
# ax[1,1].set_xlim(0,1)
# ax[0,0].set_ylim(0,1.01)
# ax[0,1].set_ylim(0,1.01)
# ax[1,0].set_ylim(0,1.01)
# ax[1,1].set_ylim(1.391,1.401)

ax2.legend(loc="lower right")
plt.subplots_adjust(left=0.12, bottom=0.15, right=0.90, top=0.94, wspace=0.20, hspace=0.20)


plt.legend()
plt.show()

# plt2.legend()
# plt2.show()

fig2.savefig('output/plots/particleCombustion/combustion.pdf')