import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
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
fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(8,11*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
fig2, ax2 = plt.subplots(nrows=4,ncols=1,figsize=(8,14*yratio))  # ,dpi=100   fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(8,8*yratio))
plt.subplots_adjust(left=0.14, bottom=0.1, right=0.90, top=0.94, wspace=0.20, hspace=0.20)
# fig = 
# ax1 = fig.add_subplot(121)
# ax2 = fig.add_subplot(122)

length = 1
condition = 2;
yIndex = 4
plot = 2    # 1 for only one plot, 2 for plot with all cases

folder = '1DflameConfined'
if length == 1:
    directory =  '900/field'
    time = 45000
    # time = 30000
elif length == 2:
    directory = 'domain768/field'
    time = 65000
    # time = 50000
elif length == 3:
    directory = 'domain1024/field'
    time = 85000
    # time = 70000

# time = 20000
if plot == 1:

    data1 = np.loadtxt('output/txt/'+folder+'/'+directory+'/'+str(time)+'.txt')
    data1 = data1[data1[:, 0].argsort()]

    position = data1[:,0]
    velocity = data1[:,4] - np.mean(data1[:,4])    # velocity signal centered at 0
    velocity = velocity/np.max(velocity)           # normalized velocity signal

    ax[0].plot(position,velocity)
    ax[0].set_ylabel('$u_\mathrm{g}\;[\mathrm{m/s}]$', fontsize=16)
    ax[0].set_xlabel('$x\;[\mathrm{m}]$', fontsize=16)

    if length == 1:
        x = np.linspace(0,0.00512,512)
    elif length == 2:
        x = np.linspace(0,0.00768,768)
    elif length == 3:
        x = np.linspace(0,0.01024,1024)


    velocityFFT = np.fft.rfft(velocity)
    print(type(velocityFFT))
    spfreq = np.fft.rfftfreq(len(velocity),x[2]-x[1])
    print(x)
    print(spfreq)
    ax[1].plot(spfreq, velocityFFT.real)
    ax[1].set_xlabel('$\mathrm{Spatial\;frequency}\;[\mathrm{1/m}]$', fontsize=16)
    ax[1].set_ylabel('$\mathrm{Amplitude}\;[\mathrm{-}]$', fontsize=16)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(-1,1))


    plt.show()

else:
    for i in range(3):
        if i == 0:
            directory =  '900/field'
            text = '0.00512'
            time = 45000
        elif i == 1:
            directory = 'domain768/field'
            text = '0.00768'
            time = 65000
        elif i == 2:
            directory = 'domain1024/field'
            text = '0.1024'
            time = 85000

        data1 = np.loadtxt('output/txt/'+folder+'/'+directory+'/'+str(time)+'.txt')
        data1 = data1[data1[:, 0].argsort()]

        position = data1[:,0]
        velocity = data1[:,4] - np.mean(data1[:,4])    # velocity signal centered at 0
        velocity = velocity/np.max(velocity)           # normalized velocity signal

        ax2[i].plot(position,velocity,color=colors[i], label='$L_x='+text+'\;\mathrm{m}$')
        ax2[i].set_ylabel('$u_\mathrm{g}\;[\mathrm{m/s}]$', fontsize=16)
        ax2[i].set_xlabel('$x\;[\mathrm{m}]$', fontsize=16)

        if i == 0:
            x = np.linspace(0,0.00512,512)
        elif i == 1:
            x = np.linspace(0,0.00768,768)
        elif i == 2:
            x = np.linspace(0,0.01024,1024)

        ax2[i].set_xlim([0,0.01024])

        velocityFFT = np.fft.rfft(velocity)
        print(type(velocityFFT))
        spfreq = np.fft.rfftfreq(len(velocity),x[2]-x[1])
        print(x)
        print(spfreq)
        ax2[3].plot(spfreq, velocityFFT.real,color=colors[i*2])
        ax2[3].set_xlabel('$\mathrm{Spatial\;frequency}\;[\mathrm{1/m}]$', fontsize=16)
        ax2[3].set_ylabel('$\mathrm{Amplitude}\;[\mathrm{-}]$', fontsize=16)
        ax2[3].set_xlim([0,2e4])
        plt.ticklabel_format(style='sci', axis='x', scilimits=(-1,1))

    ax2[0].get_xaxis().set_visible(False)
    ax2[1].get_xaxis().set_visible(False)
    ax2[2].set_xlabel('$x\;[\mathrm{m}]$', fontsize=16)
    plt.show()




# if yIndex == 1:
#     ax[0].plot(data0[256:307,0]-0.00256,data0[256:307,yIndex],c='black',linewidth=2,label='$t=t_0$') 
#     ax[0].plot(data0[308:767,0]-0.00256,data0[308:767,yIndex],c='black',linewidth=2) 
#     ax[0].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2,label='$\mathrm{AFT,Fe}$-$\mathrm{to}$-$\mathrm{FeO}$') 
#     ax[0].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2,label='$\mathrm{AFT,equilibrium}$') 

#     ax[1].plot(data1[256:767,0]-0.00256,data1[256:767,yIndex],c='black',linewidth=2,label='$t=15.0\;\mathrm{ms}$') 
#     ax[1].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2) 
#     ax[1].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2)

#     ax[2].plot(data2[256:767,0]-0.00256,data2[256:767,yIndex],c='black',linewidth=2,label='$t=30.0\mathrm{ms}$') 
#     ax[2].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2) 
#     ax[2].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2) 
    
#     ax[3].plot(data3[256:767,0]-0.00256,data3[256:767,yIndex],c='black',linewidth=2,label='$t=45.0\;\mathrm{ms}$') 
#     ax[3].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2) 
#     ax[3].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2) 
    
#     ax[4].plot(data4[256:767,0]-0.00256,data4[256:767,yIndex],c='black',linewidth=2,label='$t=60.0\;\mathrm{ms}$') 
#     ax[4].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2331,c='red',linestyle='dashed',linewidth=2) 
#     ax[4].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2) 
    
#     ax[0].set_ylim(0,2600)
#     ax[1].set_ylim(0,2600)
#     ax[2].set_ylim(0,2600)
#     ax[3].set_ylim(0,2600)
#     ax[4].set_ylim(0,2600)

#     ax[0].set_xlim(0,0.00512)
#     ax[1].set_xlim(0,0.00512)
#     ax[2].set_xlim(0,0.00512)
#     ax[3].set_xlim(0,0.00512)
#     ax[4].set_xlim(0,0.00512)

#     # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
#     # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
#     fig.supylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
#     ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
#     # fig.suptitle('Figure')
#     ax[0].legend(ncol=1, loc=(0.655, 0.125), fontsize = 14, frameon=False)
#     ax[1].legend(ncol=1, loc="lower left", fontsize = 14, frameon=False)
#     ax[2].legend(ncol=1, loc="lower left", fontsize = 14, frameon=False)
#     ax[3].legend(ncol=1, loc="lower left", fontsize = 14, frameon=False)
#     ax[4].legend(ncol=1, loc="lower left", fontsize = 14, frameon=False)

#     ax[0].get_xaxis().set_visible(False)
#     ax[1].get_xaxis().set_visible(False)
#     ax[2].get_xaxis().set_visible(False)
#     ax[3].get_xaxis().set_visible(False)

#     plt.show()

#     fig.savefig('output/plots/flame/thermalStructurePhi1.pdf')

# elif yIndex == 2:
#     ax[0].plot(data0[256:307,0]-0.00256,data0[256:307,yIndex],c='black',linewidth=3,label='$t=t_0$') 
#     ax[0].plot(data0[308:767,0]-0.00256,data0[308:767,yIndex],c='black',linewidth=3) 
#     ax[0].plot(data0[256:767,0]-0.00256,data0[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3,label='$Y_\mathrm{O_2,0}$') 

#     ax[1].plot(data1[256:767,0]-0.00256,data1[256:767,yIndex],c='black',linewidth=3,label='$t=15.0\;\mathrm{ms}$') 
#     ax[1].plot(data1[256:767,0]-0.00256,data1[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

#     ax[2].plot(data2[256:767,0]-0.00256,data2[256:767,yIndex],c='black',linewidth=3,label='$t=30.0\mathrm{ms}$') 
#     ax[2].plot(data2[256:767,0]-0.00256,data2[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

#     ax[3].plot(data3[256:767,0]-0.00256,data3[256:767,yIndex],c='black',linewidth=3,label='$t=45.0\;\mathrm{ms}$') 
#     ax[3].plot(data3[256:767,0]-0.00256,data3[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

#     ax[4].plot(data4[256:767,0]-0.00256,data4[256:767,yIndex],c='black',linewidth=3,label='$t=60.0\;\mathrm{ms}$') 
#     ax[4].plot(data4[256:767,0]-0.00256,data4[256:767,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

#     ax[0].set_ylim(0,0.232917511457580)
#     ax[1].set_ylim(0,0.232917511457580)
#     ax[2].set_ylim(0,0.232917511457580)
#     ax[3].set_ylim(0,0.232917511457580)
#     ax[4].set_ylim(0,0.232917511457580)

#     ax[0].set_xlim(0,0.00512)
#     ax[1].set_xlim(0,0.00512)
#     ax[2].set_xlim(0,0.00512)
#     ax[3].set_xlim(0,0.00512)
#     ax[4].set_xlim(0,0.00512)

#     # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
#     # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
#     fig.supylabel(r'$Y_\mathrm{O_2}\;[\mathrm{-}]$', fontsize=20)
#     ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
#     # fig.suptitle('Figure')
#     ax[0].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
#     ax[1].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
#     ax[2].legend(ncol=1, loc="upper left", fontsize = 14, frameon=False)
#     ax[3].legend(ncol=1, loc="upper left", fontsize = 14, frameon=False)
#     ax[4].legend(ncol=1, loc="upper left", fontsize = 14, frameon=False)

#     ax[0].get_xaxis().set_visible(False)
#     ax[1].get_xaxis().set_visible(False)
#     ax[2].get_xaxis().set_visible(False)
#     ax[3].get_xaxis().set_visible(False)

#     plt.show()

#     fig.savefig('output/plots/flame/O2massFracPhi1.pdf')


# add = sum(data4[456:656,1])
# avg = add/200
# print(avg)

# print(max(data1[:,1]))