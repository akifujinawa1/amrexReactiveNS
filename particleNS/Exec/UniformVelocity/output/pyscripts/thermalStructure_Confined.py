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

length = 2

if length == 1:
    directory =  'output/txt/1DflameConfined/900/field'
elif length == 2:
    directory = 'output/txt/1DflameConfined/domain768/field'
elif length == 3:
    directory = 'output/txt/1DflameConfined/domain1024/field'


# if length == 1:
#     directory =  '/../../../../mnt/D/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft512/field/'
# elif length == 2:
#     directory = '/../../../../mnt/D/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft768/field/'
# elif length == 3:
#     directory = '/../../../../mnt/D/codeBackup/flame_txt_/isochoric_Results/domainLength/fftResults/fft1024/field/'

# times = [15000,30000,45000,60000,75000]
times = [10000,20000,30000,40000,50000]
# times = [1000,3000,5000,7000,9000]
# times = [1000000,2000000,3000000,4000000,5000000]



data0 = np.loadtxt(directory+'/'+str(times[0])+'.txt')
data1 = np.loadtxt(directory+'/'+str(times[1])+'.txt')
data2 = np.loadtxt(directory+'/'+str(times[2])+'.txt')
data3 = np.loadtxt(directory+'/'+str(times[3])+'.txt')
data4 = np.loadtxt(directory+'/'+str(times[4])+'.txt')
# data5 = np.loadtxt('output/txt/1DflameConfined/'+directory+'/89900.txt')


data0 = data0[data0[:, 0].argsort()]
data1 = data1[data1[:, 0].argsort()]
data2 = data2[data2[:, 0].argsort()]
data3 = data3[data3[:, 0].argsort()]
data4 = data4[data4[:, 0].argsort()]
# data5 = data5[data5[:, 0].argsort()]

# for i in range(512):
#     data0[i,2] = 0.232917511457580
#     if i < 256:
#         data0[i,1] = 300
#     elif i < 307.2:
#         data0[i,1] = 1270
#     else:
#         data0[i,1] = 300

# plotting either temperature or mass fraction or pressure 
yIndex = 3

if yIndex == 1:
    ax[0].plot(data0[0:51,0],data0[0:51,yIndex],c='black',linewidth=2,label='$t='+str(times[0]/1e3)+'\;\mathrm{ms}$') 
    ax[0].plot(data0[52:511,0],data0[52:511,yIndex],c='black',linewidth=2) 
    ax[0].plot(data0[:,0],data0[:,yIndex]*0+2588,c='red',linestyle='dashed',linewidth=2,label='$\mathrm{AFT,Fe}$-$\mathrm{to}$-$\mathrm{FeO}$') 
    ax[0].plot(data0[:,0],data0[:,yIndex]*0+300,c=colors[0],linestyle='dotted',linewidth=2)
    # ax[0].plot(data0[:,0],data0[:,yIndex]*0+2230,c=colors[0],linestyle='dotted',linewidth=2,label='$\mathrm{AFT,equilibrium}$') 

    ax[1].plot(data1[:,0],data1[:,yIndex],c='black',linewidth=2,label='$t='+str(times[1]/1e3)+'\;\mathrm{ms}$') 
    ax[1].plot(data0[:,0],data0[:,yIndex]*0+2588,c='red',linestyle='dashed',linewidth=2) 
    ax[1].plot(data0[:,0],data0[:,yIndex]*0+300,c=colors[0],linestyle='dotted',linewidth=2)

    ax[2].plot(data2[:,0],data2[:,yIndex],c='black',linewidth=2,label='$t='+str(times[2]/1e3)+'\;\mathrm{ms}$') 
    ax[2].plot(data0[:,0],data0[:,yIndex]*0+2588,c='red',linestyle='dashed',linewidth=2) 
    ax[2].plot(data0[:,0],data0[:,yIndex]*0+300,c=colors[0],linestyle='dotted',linewidth=2) 
    
    ax[3].plot(data3[:,0],data3[:,yIndex],c='black',linewidth=2,label='$t='+str(times[3]/1e3)+'\;\mathrm{ms}$') 
    ax[3].plot(data0[:,0],data0[:,yIndex]*0+2588,c='red',linestyle='dashed',linewidth=2) 
    ax[3].plot(data0[:,0],data0[:,yIndex]*0+300,c=colors[0],linestyle='dotted',linewidth=2) 
    
    ax[4].plot(data4[:,0],data4[:,yIndex],c='black',linewidth=2,label='$t='+str(times[4]/1e3)+'\;\mathrm{ms}$') 
    ax[4].plot(data0[:,0],data0[:,yIndex]*0+2588,c='red',linestyle='dashed',linewidth=2) 
    ax[4].plot(data0[:,0],data0[:,yIndex]*0+300,c=colors[0],linestyle='dotted',linewidth=2) 
    
    ax[0].set_ylim(0,2600)
    ax[1].set_ylim(0,2600)
    ax[2].set_ylim(0,2600)
    ax[3].set_ylim(0,2600)
    ax[4].set_ylim(0,2600)
    

    # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
    # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
    fig.supylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
    ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
    # fig.suptitle('Figure')
    ax[0].legend(ncol=1, loc=(0.655, 0.125), fontsize = 14, frameon=False)
    ax[1].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[2].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[3].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[4].legend(ncol=1, loc="best", fontsize = 14, frameon=False)

    ax[0].get_xaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    ax[2].get_xaxis().set_visible(False)
    ax[3].get_xaxis().set_visible(False)


    # fig.savefig('output/plots/flame/thermalStructurePhi1.pdf')

elif yIndex == 2:
    ax[0].plot(data0[0:51,0],data0[0:51,yIndex],c='black',linewidth=3,label='$t='+str(times[0]/1e3)+'\;\mathrm{ms}$') 
    ax[0].plot(data0[52:511,0],data0[52:511,yIndex],c='black',linewidth=3) 
    ax[0].plot(data0[:,0],data0[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3,label='$Y_\mathrm{O_2,0}$') 

    ax[1].plot(data1[:,0],data1[:,yIndex],c='black',linewidth=3,label='$t='+str(times[1]/1e3)+'\;\mathrm{ms}$') 
    ax[1].plot(data1[:,0],data1[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[2].plot(data2[:,0],data2[:,yIndex],c='black',linewidth=3,label='$t='+str(times[2]/1e3)+'\;\mathrm{ms}$') 
    ax[2].plot(data2[:,0],data2[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[3].plot(data3[:,0],data3[:,yIndex],c='black',linewidth=3,label='$t='+str(times[3]/1e3)+'\;\mathrm{ms}$') 
    ax[3].plot(data3[:,0],data3[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[4].plot(data4[:,0],data4[:,yIndex],c='black',linewidth=3,label='$t='+str(times[4]/1e3)+'\;\mathrm{ms}$') 
    ax[4].plot(data4[:,0],data4[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[0].set_ylim(0,0.232917511457580)
    ax[1].set_ylim(0,0.232917511457580)
    ax[2].set_ylim(0,0.232917511457580)
    ax[3].set_ylim(0,0.232917511457580)
    ax[4].set_ylim(0,0.232917511457580)



    # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
    # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
    fig.supylabel(r'$Y_\mathrm{O_2}\;[\mathrm{-}]$', fontsize=20)
    ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
    # fig.suptitle('Figure')
    ax[0].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[1].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[2].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[3].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[4].legend(ncol=1, loc="best", fontsize = 14, frameon=False)

    ax[0].get_xaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    ax[2].get_xaxis().set_visible(False)
    ax[3].get_xaxis().set_visible(False)


    # fig.savefig('output/plots/flame/O2massFracPhi1.pdf')
elif yIndex == 3:
    ax[0].plot(data0[0:51,0],data0[0:51,yIndex],c='black',linewidth=3,label='$t='+str(times[0]/1e3)+'\;\mathrm{ms}$') 
    ax[0].plot(data0[52:511,0],data0[52:511,yIndex],c='black',linewidth=3) 
    # ax[0].plot(data0[:,0],data0[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3,label='$Y_\mathrm{O_2,0}$') 

    ax[1].plot(data1[:,0],data1[:,yIndex],c='black',linewidth=3,label='$t='+str(times[1]/1e3)+'\;\mathrm{ms}$') 
    # ax[1].plot(data1[:,0],data1[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[2].plot(data2[:,0],data2[:,yIndex],c='black',linewidth=3,label='$t='+str(times[2]/1e3)+'\;\mathrm{ms}$') 
    # ax[2].plot(data2[:,0],data2[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[3].plot(data3[:,0],data3[:,yIndex],c='black',linewidth=3,label='$t='+str(times[3]/1e3)+'\;\mathrm{ms}$') 
    # ax[3].plot(data3[:,0],data3[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[4].plot(data4[:,0],data4[:,yIndex],c='black',linewidth=3,label='$t='+str(times[4]/1e3)+'\;\mathrm{ms}$') 
    # ax[4].plot(data4[:,0],data4[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    # ax[0].set_ylim(0,6e5)
    # ax[1].set_ylim(0,6e5)
    # ax[2].set_ylim(0,6e5)
    # ax[3].set_ylim(0,6e5)
    # ax[4].set_ylim(0,6e5)

    

    # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
    # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
    fig.supylabel(r'$p\;[\mathrm{Pa}]$', fontsize=20)
    ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
    # fig.suptitle('Figure')
    ax[0].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[1].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[2].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[3].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[4].legend(ncol=1, loc="best", fontsize = 14, frameon=False)

    ax[0].get_xaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    ax[2].get_xaxis().set_visible(False)
    ax[3].get_xaxis().set_visible(False)


elif yIndex == 4:
    ax[0].plot(data0[:,0],data0[:,yIndex],c='black',linewidth=2,label='$t='+str(times[0]/1e3)+'\;\mathrm{ms}$') 
    # ax[0].plot(data0[52:511,0],data0[52:511,yIndex],c='black',linewidth=2) 
    # ax[0].plot(data0[:,0],data0[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3,label='$Y_\mathrm{O_2,0}$') 

    ax[1].plot(data1[:,0],data1[:,yIndex],c='black',linewidth=2,label='$t='+str(times[1]/1e3)+'\;\mathrm{ms}$') 
    # ax[1].plot(data1[:,0],data1[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[2].plot(data2[:,0],data2[:,yIndex],c='black',linewidth=2,label='$t='+str(times[2]/1e3)+'\;\mathrm{ms}$') 
    # ax[2].plot(data2[:,0],data2[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[3].plot(data3[:,0],data3[:,yIndex],c='black',linewidth=2,label='$t='+str(times[3]/1e3)+'\;\mathrm{ms}$') 
    # ax[3].plot(data3[:,0],data3[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[4].plot(data4[:,0],data4[:,yIndex],c='black',linewidth=2,label='$t='+str(times[4]/1e3)+'\;\mathrm{ms}$') 
    # ax[4].plot(data4[:,0],data4[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    # ax[0].set_ylim(0,6e5)
    # ax[1].set_ylim(0,6e5)
    # ax[2].set_ylim(0,6e5)
    # ax[3].set_ylim(0,6e5)
    # ax[4].set_ylim(0,6e5)


    # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
    # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
    fig.supylabel(r'$u_\mathrm{g}\;[\mathrm{m/s}]$', fontsize=20)
    ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
    # fig.suptitle('Figure')
    ax[0].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[1].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[2].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[3].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[4].legend(ncol=1, loc="best", fontsize = 14, frameon=False)

    ax[0].get_xaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    ax[2].get_xaxis().set_visible(False)
    ax[3].get_xaxis().set_visible(False)

elif yIndex == 5:
    ax[0].plot(data1[:,0],data1[:,yIndex],c='black',linewidth=2,label='$t='+str(times[0]/1e3)+'\;\mathrm{ms}$') 
    # ax[0].plot(data0[52:511,0],data0[52:511,yIndex],c='black',linewidth=2) 
    # ax[0].plot(data0[:,0],data0[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3,label='$Y_\mathrm{O_2,0}$') 

    ax[1].plot(data1[:,0],data1[:,yIndex],c='black',linewidth=2,label='$t='+str(times[1]/1e3)+'\;\mathrm{ms}$') 
    # ax[1].plot(data1[:,0],data1[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[2].plot(data2[:,0],data2[:,yIndex],c='black',linewidth=2,label='$t='+str(times[2]/1e3)+'\;\mathrm{ms}$') 
    # ax[2].plot(data2[:,0],data2[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[3].plot(data3[:,0],data3[:,yIndex],c='black',linewidth=2,label='$t='+str(times[3]/1e3)+'\;\mathrm{ms}$') 
    # ax[3].plot(data3[:,0],data3[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    ax[4].plot(data4[:,0],data4[:,yIndex],c='black',linewidth=2,label='$t='+str(times[4]/1e3)+'\;\mathrm{ms}$') 
    # ax[4].plot(data4[:,0],data4[:,yIndex]*0+0.232917511457580,c='red',linestyle='dashed',linewidth=3) 

    # ax[0].set_ylim(0,6e5)
    # ax[1].set_ylim(0,6e5)
    # ax[2].set_ylim(0,6e5)
    # ax[3].set_ylim(0,6e5)
    # ax[4].set_ylim(0,6e5)



    # ax.set_ylabel(r'$T_\mathrm{g}\;[\mathrm{K}]$', fontsize=20)
    # ax[4].set_xlabel(r'$\mathrm{x}\;[\mathrm{m}]$', fontsize=20)
    fig.supylabel(r'$\rho\;[\mathrm{kg/m^3}]$', fontsize=20)
    ax[4].set_xlabel(r'$x\;[\mathrm{m}]$', fontsize=20)
    # fig.suptitle('Figure')
    ax[0].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[1].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[2].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[3].legend(ncol=1, loc="best", fontsize = 14, frameon=False)
    ax[4].legend(ncol=1, loc="best", fontsize = 14, frameon=False)

    ax[0].get_xaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    ax[2].get_xaxis().set_visible(False)
    ax[3].get_xaxis().set_visible(False)

if length == 1:
    ax[0].set_xlim(0,0.00512)
    ax[1].set_xlim(0,0.00512)
    ax[2].set_xlim(0,0.00512)
    ax[3].set_xlim(0,0.00512)
    ax[4].set_xlim(0,0.00512)
elif length == 2:
    ax[0].set_xlim(0,0.00768)
    ax[1].set_xlim(0,0.00768)
    ax[2].set_xlim(0,0.00768)
    ax[3].set_xlim(0,0.00768)
    ax[4].set_xlim(0,0.00768)
elif length == 3:
    ax[0].set_xlim(0,0.01024)
    ax[1].set_xlim(0,0.01024)
    ax[2].set_xlim(0,0.01024)
    ax[3].set_xlim(0,0.01024)
    ax[4].set_xlim(0,0.01024)

plt.show()

# add = sum(data4[456:656,1])
# avg = add/200
# print(avg)