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

# fig,ax = plt.subplots(figsize=(10,6))
fig2,ax2 = plt.subplots(figsize=(10,6))


dp   = 27.28e-6
rp = dp/2
k  = 0.02624304274
cp = 9.739588475965419e+02
mtot = 8.363157989e-11
pi = np.pi
c = 4*pi*rp*k/(cp*mtot)

c1 = 1000-300


N = 1000000
T = 0.1
dt = T/N

time = np.arange(0,0.1,T/N)
Texact  = np.zeros(N)
# xp  = np.zeros(N)


for iT, t in enumerate(time):
    
    Texact[iT] = c1*np.exp(-c*t) + 300
    



data  = np.loadtxt('output/txt/particleDrag/data.txt')


# ax.plot(arr,y0,c="red",label='$t=0$')
# # ax.plot(arr,y1[:,0],c="green",label='$sim,t=100$')
# ax.scatter(data25[:,0],data25[:,5],c="black",marker=".",linewidths="1.5",label='$\mathrm{Numerical},t=25$')
# ax.plot(arr,y1[:,0],"--",c="blue",label='$\mathrm{Analytical},t=25$')

# # ax.scatter(data[:,0],data[:,5],c="black",label='$numerical,t=50$')
# # ax.plot(arr,y1[:,1],c="blue",label='$analytical,t=50$')

# ax.scatter(data100[:,0],data100[:,5],c="black",marker="x",linewidths="1.5",label='$\mathrm{Numerical},t=100$')
# ax.plot(arr,y1[:,1],":",c="blue",label='$\mathrm{Analytical},t=100$')

# ax.scatter(data250[:,0],data250[:,5],c="black",label='$numerical,t=250$')
# ax.plot(arr,y1[:,2],c="blue",label='$analytical,t=250$')


ax2.plot(time,Texact,c="red",label='$\mathrm{Analytical, Fixed\;}c_\mathrm{p}$')
# ax.plot(arr,y1[:,0],c="green",label='$sim,t=100$')
ax2.scatter(data[:,0],data[:,4],linewidths="0.05",c="black",label='$\mathrm{Euler\;dt=2e}$-$5,\;\mathrm{Fixed\;}c_\mathrm{p}$') # marker=".",linewidths="0.5"
ax2.scatter(data[:,0],data[:,3],linewidths="0.05",c="blue",label='$\mathrm{Euler\;dt=2e}$-$5,\;c_\mathrm{p}(T_\mathrm{p})$') # marker=".",linewidths="0.5"

# ax.scatter(data[:,0],data[:,5],c="black",label='$numerical,t=50$')
# ax.plot(arr,y1[:,1],c="blue",label='$analytical,t=50$')



# tick_spacing = 1

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

ax2.set_ylim(300,1000)

# # set axes limits
ax2.set_xlim(0,0.1)
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

fig2.savefig('output/plots/particleHeat/temperature.pdf')