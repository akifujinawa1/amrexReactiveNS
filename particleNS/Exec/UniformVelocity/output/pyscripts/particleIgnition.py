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




data1  = np.loadtxt('output/txt/particleIgnition/data1058.txt')
data2  = np.loadtxt('output/txt/particleIgnition/data1073.txt')
data3  = np.loadtxt('output/txt/particleIgnition/data1079.txt')
data4  = np.loadtxt('output/txt/particleIgnition/data1080.txt')


ax2.plot(data4[:,0],data4[:,1],c=color[0],label='$T_\mathrm{p,0}=1080\;\mathrm{K}$')
ax2.plot(data3[:,0],data3[:,1],c=color[1],label='$T_\mathrm{p,0}=1079\;\mathrm{K}$')
ax2.plot(data2[:,0],data2[:,1],c=color[2],label='$T_\mathrm{p,0}=1073\;\mathrm{K}$')
ax2.plot(data1[:,0],data1[:,1],c=color[3],label='$T_\mathrm{p,0}=1058\;\mathrm{K}$')

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

ax2.set_ylim(1050,1600)

# # set axes limits
ax2.set_xlim(0,0.03)
# ax[0,1].set_xlim(0,1)
# ax[1,0].set_xlim(0,1)
# ax[1,1].set_xlim(0,1)
# ax[0,0].set_ylim(0,1.01)
# ax[0,1].set_ylim(0,1.01)
# ax[1,0].set_ylim(0,1.01)
# ax[1,1].set_ylim(1.391,1.401)

ax2.legend(loc="lower right")


plt.legend()
plt.show()

# plt2.legend()
# plt2.show()

# fig2.savefig('output/plots/particleIgnition/ignition.pdf')