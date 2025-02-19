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

colors = ['#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D']
color = colors
fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))

data0  = np.loadtxt('output/txt/1Dflame/isobaric/field/00.txt')
data1  = np.loadtxt('output/txt/1Dflame/isobaric/field/1500.txt')
data2  = np.loadtxt('output/txt/1Dflame/isobaric/field/3000.txt')
data3  = np.loadtxt('output/txt/1Dflame/isobaric/field/4500.txt')
data4  = np.loadtxt('output/txt/1Dflame/isobaric/field/6000.txt')
data5  = np.loadtxt('output/txt/1Dflame/isobaric/field/7500.txt')
data6  = np.loadtxt('output/txt/1Dflame/isobaric/field/9000.txt')


ax2[0].scatter(data0[:,0],data0[:,1],c='black',s=5,label='$t=t_0\;\mathrm{ms}$') 
ax2[0].scatter(data1[:,0],data1[:,1],c=color[0],s=5,label='$t=1.5\mathrm{ms}$') 
ax2[0].scatter(data2[:,0],data2[:,1],c=color[1],s=5,label='$t=3.0\mathrm{ms}$') 
ax2[0].scatter(data3[:,0],data3[:,1],c=color[2],s=5,label='$t=4.5\mathrm{ms}$') 
ax2[0].scatter(data4[:,0],data4[:,1],c=color[3],s=5,label='$t=6.0\mathrm{ms}$') 
ax2[0].scatter(data5[:,0],data5[:,1],c=color[4],s=5,label='$t=7.5\mathrm{ms}$') 
ax2[0].scatter(data6[:,0],data6[:,1],c=color[5],s=5,label='$t=9.0\mathrm{ms}$') 

ax2[1].scatter(data0[:,0],data0[:,2],c='black',s=5,label='$t=t_0\;\mathrm{ms}$') 
ax2[1].scatter(data1[:,0],data1[:,2],c=color[0],s=5,label='$t=0.5\mathrm{ms}$') 
ax2[1].scatter(data2[:,0],data2[:,2],c=color[1],s=5,label='$t=1.0\mathrm{ms}$') 
ax2[1].scatter(data3[:,0],data3[:,2],c=color[2],s=5,label='$t=1.5\mathrm{ms}$') 
ax2[1].scatter(data4[:,0],data4[:,2],c=color[3],s=5,label='$t=2.0\mathrm{ms}$') 
ax2[1].scatter(data5[:,0],data5[:,2],c=color[4],s=5,label='$t=2.5\mathrm{ms}$') 
ax2[1].scatter(data6[:,0],data6[:,2],c=color[5],s=5,label='$t=3.0\mathrm{ms}$') 


ax2[0].set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
ax2[1].set_ylabel(r'$Y_\mathrm{O_2}$', fontsize=20)

ax2[1].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)

# ax2.set_ylim(1200,2750)
# ax2.set_xlim(0,0.035)

# # set axes limits
# ax2.set_xlim(0,0.01)
# ax[0,1].set_xlim(0,1)
# ax[1,0].set_xlim(0,1)
# ax[1,1].set_xlim(0,1)
# ax[0,0].set_ylim(0,1.01)
# ax[0,1].set_ylim(0,1.01)
# ax[1,0].set_ylim(0,1.01)
# ax[1,1].set_ylim(1.391,1.401)

ax2[0].legend(ncol=1, loc="top right", fontsize = 12)


# plt.legend()
plt.show()

# plt2.legend()
# plt2.show()

# fig2.savefig('output/plots/particleCombustion/combustion.pdf')