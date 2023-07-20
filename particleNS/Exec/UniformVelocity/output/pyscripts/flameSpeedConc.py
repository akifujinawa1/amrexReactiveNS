import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib

import os
import sys

from moviepy.editor import *
from moviepy.video.io.bindings import mplfig_to_npimage

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

colors = ['#000000', '#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D', '#f3d00d']
markers = ['o', 'v', '^', '<', '>', 's', 'p', '*']

# matplot subplot
fig, ax = plt.subplots(figsize=(8,8*yratio))  #fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))nrows=2,ncols=1,,dpi=100
plt.subplots_adjust(left=0.125, bottom=0.138, right=0.96, top=0.96, wspace=0.20, hspace=0.20)

time = numpy.empty((35, 8))
location = numpy.empty((35, 8))
flameSpeed = numpy.empty(8)
concArray = [600,700,800,900,1000,1100,1200,1300]

for i in range(8):
    concentration = (i+6)*100
    directory = 'output/txt/1Dflame/'+str(concentration)+'/particle/'
    Np = 0
        
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            data = np.loadtxt(f)
            # 0=time, 1=x, 2=mFe, 3=mFeO, 4=mFe3O4, 5=Tp, 6=regime
            # print(len(data[:,0]))s
            for k in range(len(data[:,0])):
                length = len(data[:,0])
                if (data[length-(k+1),6]-data[length-(k+2),6] == 1):
                    time[Np][i]=(data[length-k,0])
                    location[Np][i]=(data[length-k,1])
                    # print(time)
                    # print(location)
                    break
            Np += 1
    time[:,i].sort()
    location[:,i].sort()

    print(Np)
    # print(time[:,i])
    # print(location[:,i])

    end = len(time[:,i])
    start = len(time[:,i])-10

    A = np.vstack([time[start:end,i], np.ones(len(time[start:end,i]))]).T
    m, c = np.linalg.lstsq(A, location[start:end,i], rcond=None)[0]
    flameSpeed[i] = m*1e2
    print('Flame speed estimate for conc. = ',concentration,' is ',m*1e2,' cm/s')

    ax.scatter(time[:,i],100*location[:,i],c=colors[i],s=15,marker=markers[i],label='$'+str(concentration)+'\;\mathrm{g/cm^3}$') #s=3
    # ax.plot(time[:,i],m*time[:,i]+c,c=colors[3],linewidth=4,label=str(concentration)) #s=3


# plt.yscale("log")
# plt.xscale("log")

# ax.ticklabel_format(useOffset=False)
ax.set_ylabel(r'$x\;[\mathrm{cm}]$', fontsize=20)
ax.set_xlabel(r'$\mathrm{time}\;[\mathrm{s}]$', fontsize=20)



# ax.scatter(concArray,flameSpeed,c='black',s=3,label='$\mathrm{Flame\;speed}$')



# # # ax2.set_ylim(1200,2750)
# # # ax2.set_xlim(0,0.035)

# # # # set axes limits
# # # ax2.set_xlim(0,0.01)
# # # ax[0,1].set_xlim(0,1)
# # # ax[1,0].set_xlim(0,1)
# # # ax[1,1].set_xlim(0,1)
# # # ax[0,0].set_ylim(0,1.01)
# # # ax[0,1].set_ylim(0,1.01)
# # # ax[1,0].set_ylim(0,1.01)
# # # ax[1,1].set_ylim(1.391,1.401)

# ax.legend(ncol=1, loc="best", fontsize = 16)



plt.legend()
plt.show()

# # # plt2.legend()
# # # plt2.show()

fig.savefig('output/plots/flame/flameSpeedX-T.pdf')

# with open('output/txt/1Dflame/phi1/x-t.txt', 'w') as text_file:
#     for i in range(len(time)):
#         timeval = time[i]
#         locationval = location[i]
#         text_file.write(str(timeval)+' '+str(locationval)+'\n')

# text_file.close()