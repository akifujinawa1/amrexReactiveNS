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


x0 = 0
L  = 1
D  = 2*10**(-5)
nlim = 100
pi   = np.pi

N = 1024
arr = np.arange(1/(2*N),1+1/(2*N),1/N)
y0  = np.zeros(N)
y1  = np.matlib.zeros((N, 3))
yInt = y0

for idx, x in enumerate(arr):
    xval = -100*(x-0.5)**2
    y0val = np.exp(xval)
    y0[idx] = y0val

t = np.array([25,100,250])
for it, time in enumerate(t):                    # loop thorugh time for plotting
    for idx, x in enumerate(arr):                # loop through x
        val = np.trapz(y0,x=arr)
        for i in range(1,50):                    # loop through n for solution
            # print(y0*np.sin(i*pi*arr))
            yInt = np.multiply(y0,np.cos(i*pi*arr))
            b    = 2*np.trapz(yInt,x=arr)        # calculate coeff b_n
            val += b*np.cos(i*pi*x)*np.exp(-i*i*pi*pi*D*time) # calculate actual solution to nth accuracy  

        # print(time)
        y1[idx,it] = val

# print(y1[256,:])

data25  = np.loadtxt('output/txt/multiGasDiffusion/time025.txt')
data100 = np.loadtxt('output/txt/multiGasDiffusion/time0100.txt')


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


ax2.plot(arr,1-y0,c="red",label='$t=0$')
# ax.plot(arr,y1[:,0],c="green",label='$sim,t=100$')
ax2.scatter(data25[:,0],data25[:,6],c="black",marker=".",linewidths="1.5",label='$\mathrm{Numerical},t=25$')
ax2.plot(arr,1-y1[:,0],"--",c="blue",label='$\mathrm{Analytical},t=25$')

# ax.scatter(data[:,0],data[:,5],c="black",label='$numerical,t=50$')
# ax.plot(arr,y1[:,1],c="blue",label='$analytical,t=50$')

ax2.scatter(data100[:,0],data100[:,6],c="black",marker="x",linewidths="1.5",label='$\mathrm{Numerical},t=100$')
ax2.plot(arr,1-y1[:,1],":",c="blue",label='$\mathrm{Analytical},t=100$')





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
ax2.set_ylabel(r'$Y_\mathrm{N_2}$', fontsize=20)
ax2.set_xlabel(r'$x$', fontsize=20)

# # set axes limits
# ax[0,0].set_xlim(0,1)
# ax[0,1].set_xlim(0,1)
# ax[1,0].set_xlim(0,1)
# ax[1,1].set_xlim(0,1)
# ax[0,0].set_ylim(0,1.01)
# ax[0,1].set_ylim(0,1.01)
# ax[1,0].set_ylim(0,1.01)
# ax[1,1].set_ylim(1.391,1.401)


plt.legend()
plt.show()

# plt2.legend()
# plt2.show()

fig2.savefig('output/plots/multiGasDiffusion/YN2.pdf')
