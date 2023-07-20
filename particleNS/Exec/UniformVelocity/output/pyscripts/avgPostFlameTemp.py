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

colors = ['#9B0909', '#7E301E', '#B42E0F', '#FE3C0F', '#fe770f', '#F35D0D']
color = colors

# tfinal = 40000

# frames = 400
# fpsval = 20
# duration = frames/fpsval
# timestep = 1/duration
# framestep = tfinal/frames
# scale = framestep/timestep
 
# matplot subplot
fig, ax = plt.subplots(nrows=5,ncols=1,figsize=(8,16*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
plt.subplots_adjust(left=0.14, bottom=0.15, right=0.90, top=0.94, wspace=0.20, hspace=0.20)
# fig = 
# ax1 = fig.add_subplot(121)
# ax2 = fig.add_subplot(122)

for i in range(8):
    concentration = (i+6)*100
    'output/txt/1Dflame/'+str(concentration)+'/field/'

    data0 = np.loadtxt('output/txt/1Dflame/'+str(concentration)+'/field/00.txt')
    data1 = np.loadtxt('output/txt/1Dflame/'+str(concentration)+'/field/7500.txt')
    data2 = np.loadtxt('output/txt/1Dflame/'+str(concentration)+'/field/15000.txt')
    data3 = np.loadtxt('output/txt/1Dflame/'+str(concentration)+'/field/22500.txt')
    data4 = np.loadtxt('output/txt/1Dflame/'+str(concentration)+'/field/30000.txt')
    data5 = np.loadtxt('output/txt/1Dflame/'+str(concentration)+'/field/37500.txt')

    data0 = data0[data0[:, 0].argsort()]
    data1 = data1[data1[:, 0].argsort()]
    data2 = data2[data2[:, 0].argsort()]
    data3 = data3[data3[:, 0].argsort()]
    data4 = data4[data4[:, 0].argsort()]
    data5 = data5[data5[:, 0].argsort()]

    for i in range(768):
        data0[i,2] = 0.232917511457580
        if i < 256:
            data0[i,1] = 300
        elif i < 307.2:
            data0[i,1] = 1270
        else:
            data0[i,1] = 300

    # plotting either temperature or mass fraction
    yIndex = 1

    add = sum(data4[456:656,1])
    avg = add/200
    print(avg)