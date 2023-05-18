import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker

import os
import sys

# sys.path.append(os.path.join('pyblish','plots'))
# import publish 


plt.close('all')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams["mathtext.fontset"] = 'cm'
plt.rc('font', family='serif', size='12')

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False


fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(10,6))
fig.tight_layout(pad=2.4)

labels = ["x","rho","v","p","eps","yO2","yN2","Tg","gamma"]

data = np.loadtxt('output/txt/multiGasSod/time0.txt')
exact = np.loadtxt('output/txt/multiGasSod/exact.txt')
ideal = np.loadtxt('output/txt/multiGasSod/idealgas.txt')

gamma = [1.4] * len(exact)


ax[0,0].scatter(ideal[:,0],ideal[:,1],marker=".",c="red",linewidths="0.1",label='$\gamma=1.4$')
ax[0,1].scatter(ideal[:,0],ideal[:,2],marker=".",c="red",linewidths="0.1",label='$\gamma=1.4$')
ax[1,0].scatter(ideal[:,0],ideal[:,3],marker=".",c="red",linewidths="0.1",label='$\gamma=1.4$')


ax[0,0].scatter(data[:,0],data[:,1],marker=".",c="blue",linewidths="0.1",label='$\mathrm{multi\;species}$')
ax[0,1].scatter(data[:,0],data[:,2],marker=".",c="blue",linewidths="0.1",label='$\mathrm{multi\;species}$')
ax[1,0].scatter(data[:,0],data[:,3],marker=".",c="blue",linewidths="0.1",label='$\mathrm{multi\;species}$')
ax[1,1].scatter(data[:,0],data[:,8],marker=".",c="blue",linewidths="0.1",label='$\mathrm{multi\;species}$')


ax[0,0].plot(exact[:,0],exact[:,1],c="black",label='$\mathrm{exact}$')
ax[0,1].plot(exact[:,0],exact[:,2],c="black",label='$\mathrm{exact}$')
ax[1,0].plot(exact[:,0],exact[:,3],c="black",label='$\mathrm{exact}$')
ax[1,1].plot(exact[:,0],gamma,c="black",label='$\mathrm{exact}$')


ax[0,0].legend(loc="upper right")
ax[0,1].legend(loc=(0.4,0.1))
ax[1,0].legend(loc="upper right")
ax[1,1].legend(loc="lower left")



tick_spacing = 1

# set labels
ax[1,0].set_xlabel(r'$x$',fontsize = 20)
ax[1,1].set_xlabel(r'$x$',fontsize = 20)
ax[0,0].set_ylabel(r'$\rho$',fontsize = 20)
ax[0,1].set_ylabel(r'$u$',fontsize = 20)
ax[1,0].set_ylabel(r'$p$',fontsize = 20)
ax[1,1].set_ylabel(r'$\gamma$',fontsize = 20)

# set axes limits
ax[0,0].set_xlim(0,1)
ax[0,1].set_xlim(0,1)
ax[1,0].set_xlim(0,1)
ax[1,1].set_xlim(0,1)
ax[0,0].set_ylim(0,1.01)
ax[0,1].set_ylim(0,1.01)
ax[1,0].set_ylim(0,1.01)
ax[1,1].set_ylim(1.391,1.401)


plt.legend()
plt.show()

fig.savefig('output/plots/multiGasEuler/sod.pdf')
