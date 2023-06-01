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


rho  = 1.171957919
rhop = 7867.533139
dp   = 2*1.364e-5
mu   = 1.843725182e-5
ug   = 5
up0  = 0
urel = ug-up0
Re   = rho*dp*urel/mu
CD   = (24/Re)*(1+0.15*Re**(0.687))
alpha  = 3*CD*rho*0.25/(dp*rhop)


N = 1000000
T = 0.1
dt = T/N

time = np.arange(0,0.1,T/N)
up  = np.zeros(N)
xp  = np.zeros(N)


for iT, t in enumerate(time):
    
    up[iT] = (-up0 -alpha*ug*ug*t + alpha*up0*ug*t)/(-alpha*ug*t+alpha*up0*t-1)
    # urel = ug-up[iT]
    # Re   = rho*dp*urel/mu
    # CD   = (24/Re)*(1+0.15*Re**(0.687))
    



data  = np.loadtxt('output/txt/particleDrag/data_fixeduRel.txt')
data2  = np.loadtxt('output/txt/particleDrag/data_dynamicuRel.txt')


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


ax2.plot(time,up,c="red",label='$\mathrm{Analytical,\;Fixed\;}u_\mathrm{rel}$')
# ax.plot(arr,y1[:,0],c="green",label='$sim,t=100$')
ax2.scatter(data[:,0],data[:,2],linewidths="0.05",c="black",label='$\mathrm{Euler\;dt=2e}$-$5,\;\mathrm{Fixed\;}u_\mathrm{rel}$') # marker=".",linewidths="0.5"
ax2.scatter(data2[:,0],data2[:,2],linewidths="0.05",c="blue",label='$\mathrm{Euler\;dt=2e}$-$5,\;u_\mathrm{rel}(t)$') # marker=".",linewidths="0.5"

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
ax2.set_ylabel(r'$u_\mathrm{p}\;[\mathrm{m/s}]$', fontsize=20)
ax2.set_xlabel(r'$\mathrm{time}\;[\mathrm{s}]$', fontsize=20)

ax2.set_ylim(0,5)

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

fig2.savefig('output/plots/particleDrag/up_time.pdf')