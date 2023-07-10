

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker

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
plt.rc('font', family='serif', size='12')

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False


fig,ax = plt.subplots(figsize=(8,8*yratio))
fig.tight_layout(pad=2.4)

labels = ["x","rho","v","p","eps","yO2","yN2","Tg","gamma"]

data64 = np.loadtxt('output/txt/multiGasSod/data64.txt')
data128 = np.loadtxt('output/txt/multiGasSod/data128.txt')
data256 = np.loadtxt('output/txt/multiGasSod/data256.txt')
data512 = np.loadtxt('output/txt/multiGasSod/data512.txt')
data1024 = np.loadtxt('output/txt/multiGasSod/data1024.txt')
data2048 = np.loadtxt('output/txt/multiGasSod/data2048.txt')

data64 = data64[data64[:, 0].argsort()]
data128 = data128[data128[:, 0].argsort()]
data256 = data256[data256[:, 0].argsort()]
data512 = data512[data512[:, 0].argsort()]
data1024 = data1024[data1024[:, 0].argsort()]
data2048 = data2048[data2048[:, 0].argsort()]

data64e = np.loadtxt('output/txt/multiGasSod/field1MUSCL64noAMR1.txt')
data128e = np.loadtxt('output/txt/multiGasSod/field1MUSCL128noAMR1.txt')
data256e = np.loadtxt('output/txt/multiGasSod/field1MUSCL256noAMR1.txt')
data512e = np.loadtxt('output/txt/multiGasSod/field1MUSCL512noAMR1.txt')
data1024e = np.loadtxt('output/txt/multiGasSod/field1MUSCL1024noAMR1.txt')
data2048e = np.loadtxt('output/txt/multiGasSod/field1MUSCL2048noAMR1.txt')

data64e = data64e[data64e[:, 0].argsort()]
data128e = data128e[data128e[:, 0].argsort()]
data256e = data256e[data256e[:, 0].argsort()]
data512e = data512e[data512e[:, 0].argsort()]
data1024e = data1024e[data1024e[:, 0].argsort()]
data2048e = data2048e[data2048e[:, 0].argsort()]

exact = np.loadtxt('output/txt/multiGasSod/exact.txt')
ideal = np.loadtxt('output/txt/multiGasSod/idealgas.txt')

gamma = [1.4] * len(exact)

iter = [0]*6
error = [0]*6
errorEuler = [0]*6
nx = [64, 128, 256, 512, 1024, 2048]

print(len(exact[:,0]))

for i in range(len(exact[:,1])):
    for j in range(6):
    # j = 0
        if (i == (iter[j]+1)*len(exact[:,0])/(2**(6+j))):
            print(i)
            avg = 0.5*(exact[i,1]+exact[i-1,1])
            if (j == 0):
                error[j] += abs(avg-data64[iter[j],1])/(2**(6+j))
                errorEuler[j] += abs(avg-data64e[iter[j],1])/(2**(6+j))
            elif (j == 1):
                error[j] += abs(avg-data128[iter[j],1])/(2**(6+j))
                errorEuler[j] += abs(avg-data128e[iter[j],1])/(2**(6+j))
            elif (j == 2):
                error[j] += abs(avg-data256[iter[j],1])/(2**(6+j))
                errorEuler[j] += abs(avg-data256e[iter[j],1])/(2**(6+j))
            elif (j == 3):
                error[j] += abs(avg-data512[iter[j],1])/(2**(6+j))
                errorEuler[j] += abs(avg-data512e[iter[j],1])/(2**(6+j))
            elif (j == 4):
                error[j] += abs(avg-data1024[iter[j],1])/(2**(6+j))
                errorEuler[j] += abs(avg-data1024e[iter[j],1])/(2**(6+j))
            else:
                error[j] += abs(avg-data2048[iter[j],1])/(2**(6+j))
                errorEuler[j] += abs(avg-data2048e[iter[j],1])/(2**(6+j))
            iter[j] += 1


print(iter[j])

print('error in 64 cell: ',error[0])
print('error in 128 cell: ',error[1])
print('error in 256 cell: ',error[2])
print('error in 512 cell: ',error[3])
print('error in 1024 cell: ',error[4])
print('error in 2048 cell: ',error[5])

first = [error[0], error[0]/2, error[0]/4, error[0]/8, error[0]/16, error[0]/32]
second = [error[0], error[0]/4, error[0]/16, error[0]/64, error[0]/256, error[0]/1024]

ax.scatter(nx,error,marker=".",c="red",linewidths="2",label='$L_1,\mathrm{O_2}$-$\mathrm{N_2}$')
ax.scatter(nx,errorEuler,marker=".",c="orange",linewidths="2",label='$L_1,\gamma =1.4$')

ax.plot(nx,first,c="black",linewidth="1.5",label='$\mathrm{O}(1)$')
ax.plot(nx,second,c="blue",linewidth="1.5",label='$\mathrm{O}(2)$')
ax.set_yscale('log')
ax.set_xscale('log')
# ax[0,0].scatter(ideal[:,0],ideal[:,1],marker=".",c="red",linewidths="0.1",label='$\gamma=1.4$')
# ax[0,1].scatter(ideal[:,0],ideal[:,2],marker=".",c="red",linewidths="0.1",label='$\gamma=1.4$')
# ax[1,0].scatter(ideal[:,0],ideal[:,3],marker=".",c="red",linewidths="0.1",label='$\gamma=1.4$')


# ax[0,0].scatter(data[:,0],data[:,1],marker=".",c="blue",linewidths="0.1",label='$\mathrm{multi\;species}$')
# ax[0,1].scatter(data[:,0],data[:,2],marker=".",c="blue",linewidths="0.1",label='$\mathrm{multi\;species}$')
# ax[1,0].scatter(data[:,0],data[:,3],marker=".",c="blue",linewidths="0.1",label='$\mathrm{multi\;species}$')
# ax[1,1].scatter(data[:,0],data[:,8],marker=".",c="blue",linewidths="0.1",label='$\mathrm{multi\;species}$')


# ax[0,0].plot(exact[:,0],exact[:,1],c="black",label='$\mathrm{exact}$')
# ax[0,1].plot(exact[:,0],exact[:,2],c="black",label='$\mathrm{exact}$')
# ax[1,0].plot(exact[:,0],exact[:,3],c="black",label='$\mathrm{exact}$')
# ax[1,1].plot(exact[:,0],gamma,c="black",label='$\mathrm{exact}$')


# ax[0,0].legend(loc="upper right")
# ax[0,1].legend(loc=(0.4,0.1))
# ax[1,0].legend(loc="upper right")

ax.legend(fontsize="20",loc="best")



# tick_spacing = 1

# # set labels
ax.set_xlabel(r'$N_x$',fontsize = 20)
# ax[1,1].set_xlabel(r'$x$',fontsize = 20)
ax.set_ylabel(r'$\epsilon$',fontsize = 20)
# ax[0,1].set_ylabel(r'$u$',fontsize = 20)
# ax[1,0].set_ylabel(r'$p$',fontsize = 20)
# ax[1,1].set_ylabel(r'$\gamma$',fontsize = 20)

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

fig.savefig('output/plots/multiGasEuler/convEuler.pdf')
