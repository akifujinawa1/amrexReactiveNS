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

colors = ['#7E301E', '#B42E0F', '#FE3C0F', '#fe770f']
color = colors

fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))

data0  = np.loadtxt('output/txt/1Dflame/isobaric/0.000100iteration10800.txt')
data1  = np.loadtxt('output/txt/1Dflame/isobaric/0.000200iteration23200.txt')
data2  = np.loadtxt('output/txt/1Dflame/isobaric/0.000300iteration36100.txt')
data3  = np.loadtxt('output/txt/1Dflame/isobaric/0.000400iteration49300.txt')
data4  = np.loadtxt('output/txt/1Dflame/isobaric/0.000500iteration62300.txt')
data5  = np.loadtxt('output/txt/1Dflame/isobaric/0.000600iteration75300.txt')
data6  = np.loadtxt('output/txt/1Dflame/isobaric/0.000700iteration88600.txt')
data7  = np.loadtxt('output/txt/1Dflame/isobaric/0.000800iteration102100.txt')
data8  = np.loadtxt('output/txt/1Dflame/isobaric/0.000900iteration115700.txt')
data9  = np.loadtxt('output/txt/1Dflame/isobaric/0.001000iteration129400.txt')

data0[:,0], data0[:,1], data0[:,2] = zip(*sorted(zip(data0[:,0], data0[:,1], data0[:,2])))
data1[:,0], data1[:,1], data1[:,2] = zip(*sorted(zip(data1[:,0], data1[:,1], data1[:,2])))
data2[:,0], data2[:,1], data2[:,2] = zip(*sorted(zip(data2[:,0], data2[:,1], data2[:,2])))
data3[:,0], data3[:,1], data3[:,2] = zip(*sorted(zip(data3[:,0], data3[:,1], data3[:,2])))
data4[:,0], data4[:,1], data4[:,2] = zip(*sorted(zip(data4[:,0], data4[:,1], data4[:,2])))
data5[:,0], data5[:,1], data5[:,2] = zip(*sorted(zip(data5[:,0], data5[:,1], data5[:,2])))
data6[:,0], data6[:,1], data6[:,2] = zip(*sorted(zip(data6[:,0], data6[:,1], data6[:,2])))
data7[:,0], data7[:,1], data7[:,2] = zip(*sorted(zip(data7[:,0], data7[:,1], data7[:,2])))
data8[:,0], data8[:,1], data8[:,2] = zip(*sorted(zip(data8[:,0], data8[:,1], data8[:,2])))
data9[:,0], data9[:,1], data9[:,2] = zip(*sorted(zip(data9[:,0], data9[:,1], data9[:,2])))


ax2[0].plot(data0[:,0],data0[:,1],c='black',linewidth=1,label='$t=0.1\mathrm{ms}$') 
ax2[0].plot(data1[:,0],data1[:,1],c=color[0],linewidth=1,label='$t=0.2\mathrm{ms}$') 
ax2[0].plot(data2[:,0],data2[:,1],c=color[1],linewidth=1,label='$t=0.3\mathrm{ms}$') 
ax2[0].plot(data3[:,0],data3[:,1],c=color[2],linewidth=1,label='$t=0.4\mathrm{ms}$') 
ax2[0].plot(data4[:,0],data4[:,1],c=color[3],linewidth=1,label='$t=0.5\mathrm{ms}$') 
ax2[0].plot(data5[:,0],data5[:,1],c='black',linewidth=1,label='$t=0.6\mathrm{ms}$') 
ax2[0].plot(data6[:,0],data6[:,1],c=color[0],linewidth=1,label='$t=0.7\mathrm{ms}$') 
ax2[0].plot(data7[:,0],data7[:,1],c=color[1],linewidth=1,label='$t=0.8\mathrm{ms}$') 
ax2[0].plot(data8[:,0],data8[:,1],c=color[2],linewidth=1,label='$t=0.9\mathrm{ms}$') 
ax2[0].plot(data9[:,0],data9[:,1],c=color[3],linewidth=1,label='$t=1.0\mathrm{ms}$') 

ax2[1].plot(data0[:,0],data0[:,2],c='black',linewidth=1,label='$t=0.1\mathrm{ms}$') 
ax2[1].plot(data1[:,0],data1[:,2],c=color[0],linewidth=1,label='$t=0.2\mathrm{ms}$') 
ax2[1].plot(data2[:,0],data2[:,2],c=color[1],linewidth=1,label='$t=0.3\mathrm{ms}$') 
ax2[1].plot(data3[:,0],data3[:,2],c=color[2],linewidth=1,label='$t=0.4\mathrm{ms}$') 
ax2[1].plot(data4[:,0],data4[:,2],c=color[3],linewidth=1,label='$t=0.5\mathrm{ms}$') 
ax2[1].plot(data5[:,0],data5[:,2],c='black',linewidth=1,label='$t=0.6\mathrm{ms}$') 
ax2[1].plot(data6[:,0],data6[:,2],c=color[0],linewidth=1,label='$t=0.7\mathrm{ms}$') 
ax2[1].plot(data7[:,0],data7[:,2],c=color[1],linewidth=1,label='$t=0.8\mathrm{ms}$') 
ax2[1].plot(data8[:,0],data8[:,2],c=color[2],linewidth=1,label='$t=0.9\mathrm{ms}$') 
ax2[1].plot(data9[:,0],data9[:,2],c=color[3],linewidth=1,label='$t=1.0\mathrm{ms}$') 


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

ax2[1].legend(loc="lower right")


# plt.legend()
plt.show()

# plt2.legend()
# plt2.show()

# fig2.savefig('output/plots/particleCombustion/combustion.pdf')