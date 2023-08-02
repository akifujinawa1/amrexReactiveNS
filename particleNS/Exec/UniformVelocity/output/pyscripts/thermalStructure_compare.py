import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import numpy.matlib

import os
import sys

# from moviepy.editor import *
# from moviepy.video.io.bindings import mplfig_to_npimage

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

colors = ['#000000', '#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D', '#f3d00d', '#9b9765']
markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'P']


# tfinal = 40000

# frames = 400
# fpsval = 20
# duration = frames/fpsval
# timestep = 1/duration
# framestep = tfinal/frames
# scale = framestep/timestep
 
# matplot subplot
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,8*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
plt.subplots_adjust(left=0.14, bottom=0.15, right=0.90, top=0.94, wspace=0.20, hspace=0.20)
# fig = 
# ax1 = fig.add_subplot(121)
# ax2 = fig.add_subplot(122)

condition = 1;
yIndex = 1

if condition == 1:
    folder = '1Dflame'
else:
    folder = '1DflameConfined'

if condition == 1:
    concentrations = [725,750,800,900,1000,1100,1200,1300,1400]
elif condition == 2:
    concentrations = [600,700,800,900,1000,1100,1200,1300,1400]

Nparams = len(concentrations)    
time = 59000

for i in range(Nparams):
    conc = concentrations[i]
    data1 = np.loadtxt('output/txt/'+folder+'/'+str(conc)+'/field/'+str(time)+'.txt')
    data1 = data1[data1[:, 0].argsort()]
    ax.plot(data1[256:767,0]-0.00256,data1[256:767,1],c=colors[i],linewidth=2,label='$'+str(conc)+'\;\mathrm{g/m^3}$') 

ax.set_ylabel('$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=16)
ax.set_xlabel('$x\;[\mathrm{m}]$', fontsize=16)
ax.legend(ncol=2, loc="best", fontsize = 16, frameon = False )

plt.show()