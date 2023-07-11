

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

data64 = np.loadtxt('output/txt/multiGasDiffusion/data64.txt')
data128 = np.loadtxt('output/txt/multiGasDiffusion/data128.txt')
data256 = np.loadtxt('output/txt/multiGasDiffusion/data256.txt')
data512 = np.loadtxt('output/txt/multiGasDiffusion/data512.txt')
data1024 = np.loadtxt('output/txt/multiGasDiffusion/data1024.txt')
data4096 = np.loadtxt('output/txt/multiGasDiffusion/data4096.txt')

data64 = data64[data64[:, 0].argsort()]
data128 = data128[data128[:, 0].argsort()]
data256 = data256[data256[:, 0].argsort()]
data512 = data512[data512[:, 0].argsort()]
data1024 = data1024[data1024[:, 0].argsort()]
data4096 = data4096[data4096[:, 0].argsort()]

# calculate exact solution here

# print(data64[:,5])

x0 = 0
L  = 0.1024
D  = 2*10**(-5)
nlim = 100
pi   = np.pi




N = 4096
arr = np.arange(L/(2*N),L+L/(2*N),L/N)
# arr = data64[:,0]
y0  = np.zeros(N)
y4096  = np.zeros(N)
yInt = y0

# print(arr)

for idx, x in enumerate(arr):
    y0val = np.exp(-5000*(x-0.0512)*(x-0.0512))
    y0[idx] = y0val

# print(y0)

t = np.array([0.0001])
for it, time in enumerate(t):                    # loop thorugh time for plotting
    for idx, x in enumerate(arr):                # loop through x
        val = np.trapz(y0,x=arr)/L               # a0 values
        for i in range(1,50):                    # loop through n for solution
            # print(y0*np.sin(i*pi*arr))
            yInt = np.multiply(y0,np.cos(i*pi*arr/L))              # integrand
            b    = 2*np.trapz(yInt,x=arr)/L                        # a_n
            val += b*np.cos(i*pi*x/L)*np.exp(-i*i*pi*pi*D*time/L**2) # calculate actual solution to nth accuracy  

        # print(time)
        y4096[idx] = val

# print(data64[:,5]-y64)
ax.scatter(data64[:,0],data64[:,5],marker=".",c="red",linewidths="0.1",label='$num.$')
ax.scatter(arr,y4096,marker=".",c="black",linewidths="0.1",label='$exact$')
plt.show()

# # print(arr-data64[:,0])

# N = 128
# # arr = np.arange(L/(2*N),L+L/(2*N),L/N)
# arr = data128[:,0]

# y0  = np.zeros(N)
# y128  = np.zeros(N)
# yInt = y0

# for idx, x in enumerate(arr):
#     # xval = -100*(x-0.5)**2
#     y0val = np.exp(-5000*(x-0.0512)*(x-0.0512))
#     y0[idx] = y0val

# t = np.array([0.1])
# for it, time in enumerate(t):                    # loop thorugh time for plotting
#     for idx, x in enumerate(arr):                # loop through x
#         val = np.trapz(y0,x=arr)/L               # a0 values
#         for i in range(1,50):                    # loop through n for solution
#             # print(y0*np.sin(i*pi*arr))
#             yInt = np.multiply(y0,np.cos(i*pi*arr/L))              # integrand
#             b    = 2*np.trapz(yInt,x=arr)/L                        # a_n
#             val += b*np.cos(i*pi*x/L)*np.exp(-i*i*pi*pi*D*time/L**2) # calculate actual solution to nth accuracy  
#         # print(time)
#         y128[idx] = val

# N = 256
# # arr = np.arange(L/(2*N),L+L/(2*N),L/N)
# arr = data256[:,0]
# y0  = np.zeros(N)
# y256  = np.zeros(N)
# yInt = y0

# for idx, x in enumerate(arr):
#     # xval = -100*(x-0.5)**2
#     y0val = np.exp(-5000*(x-0.0512)*(x-0.0512))
#     y0[idx] = y0val

# t = np.array([0.1])
# for it, time in enumerate(t):                    # loop thorugh time for plotting
#     for idx, x in enumerate(arr):                # loop through x
#         val = np.trapz(y0,x=arr)/L               # a0 values
#         for i in range(1,50):                    # loop through n for solution
#             # print(y0*np.sin(i*pi*arr))
#             yInt = np.multiply(y0,np.cos(i*pi*arr/L))              # integrand
#             b    = 2*np.trapz(yInt,x=arr)/L                        # a_n
#             val += b*np.cos(i*pi*x/L)*np.exp(-i*i*pi*pi*D*time/L**2) # calculate actual solution to nth accuracy  
#         # print(time)
#         y256[idx] = val




# N = 512
# # arr = np.arange(L/(2*N),L+L/(2*N),L/N)
# arr = data512[:,0]
# y0  = np.zeros(N)
# y512  = np.zeros(N)
# yInt = y0

# for idx, x in enumerate(arr):
#     # xval = -100*(x-0.5)**2
#     y0val = np.exp(-5000*(x-0.0512)*(x-0.0512))
#     y0[idx] = y0val

# t = np.array([0.1])
# for it, time in enumerate(t):                    # loop thorugh time for plotting
#     for idx, x in enumerate(arr):                # loop through x
#         val = np.trapz(y0,x=arr)/L               # a0 values
#         for i in range(1,50):                    # loop through n for solution
#             # print(y0*np.sin(i*pi*arr))
#             yInt = np.multiply(y0,np.cos(i*pi*arr/L))              # integrand
#             b    = 2*np.trapz(yInt,x=arr)/L                        # a_n
#             val += b*np.cos(i*pi*x/L)*np.exp(-i*i*pi*pi*D*time/L**2) # calculate actual solution to nth accuracy  

#         # print(time)
#         y512[idx] = val

# N = 1024
# # arr = np.arange(L/(2*N),L+L/(2*N),L/N)
# arr = data1024[:,0]
# y0  = np.zeros(N)
# y1024 = np.zeros(N)
# yInt = y0

# for idx, x in enumerate(arr):
#     # xval = -100*(x-0.5)**2
#     y0val = np.exp(-5000*(x-0.0512)*(x-0.0512))
#     y0[idx] = y0val

# t = np.array([0.1])
# for it, time in enumerate(t):                    # loop thorugh time for plotting
#     for idx, x in enumerate(arr):                # loop through x
#         val = np.trapz(y0,x=arr)/L               # a0 values
#         for i in range(1,50):                    # loop through n for solution
#             # print(y0*np.sin(i*pi*arr))
#             yInt = np.multiply(y0,np.cos(i*pi*arr/L))              # integrand
#             b    = 2*np.trapz(yInt,x=arr)/L                        # a_n
#             val += b*np.cos(i*pi*x/L)*np.exp(-i*i*pi*pi*D*time/L**2) # calculate actual solution to nth accuracy  

#         # print(time)
#         y1024[idx] = val

# ax.scatter(arr,y1024,marker=".",c="red",linewidths="0.1",label='$exact$')
# ax.scatter(data1024[:,0],data1024[:,5],marker=".",c="black",linewidths="0.1",label='$L_1,\mathrm{O_2}$-$\mathrm{N_2}$')
# plt.show()
# print(y128)

iter = [0]*5
error = [0]*5
nx = [64, 128, 256, 512, 1024]

# print(len(y1))

# data4096=y4096

for i in range(len(y4096)):
    for j in range(5):
    # j = 0
        if (i == (iter[j]+1)*len(y4096)/(2**(6+j))):
            # print(i)
            avg = 0.5*(y4096[i]+y4096[i-1])
            if (j == 0):
                error[j] += abs(avg-data64[iter[j],5])/(2**(6+j))
            elif (j == 1):
                error[j] += abs(avg-data128[iter[j],5])/(2**(6+j))
            elif (j == 2):
                error[j] += abs(avg-data256[iter[j],5])/(2**(6+j))
            elif (j == 3):
                error[j] += abs(avg-data512[iter[j],5])/(2**(6+j))
            else:
                error[j] += abs(avg-data1024[iter[j],5])/(2**(6+j))
            iter[j] += 1


    # for j in range(5):
    # # j = 0
    #     if (i == (iter[j]+1)*len(y1)/(2**(6+j))):
    #         # print(i)
    #         avg = 0.5*(y1[i]+y1[i-1])
    #         if (j == 0):
    #             error[j] += abs(y64-data64[iter[j],5])/(2**(6+j))
    #         elif (j == 1):
    #             error[j] += abs(avg-data128[iter[j],5])/(2**(6+j))
    #         elif (j == 2):
    #             error[j] += abs(avg-data256[iter[j],5])/(2**(6+j))
    #         elif (j == 3):
    #             error[j] += abs(avg-data512[iter[j],5])/(2**(6+j))
    #         else:
    #             error[j] += abs(avg-data1024[iter[j],5])/(2**(6+j))
    #         iter[j] += 1


# print(data64[:5])
# print(y64)

print('error in 64 cell: ',error[0])
print('error in 128 cell: ',error[1])
print('error in 256 cell: ',error[2])
print('error in 512 cell: ',error[3])
print('error in 1024 cell: ',error[4])

first = [error[0], error[0]/2, error[0]/4, error[0]/8, error[0]/16]
second = [error[0], error[0]/4, error[0]/16, error[0]/64, error[0]/256]





ax.scatter(nx,error,marker=".",c="red",linewidths="2",label='$L_1,\mathrm{O_2}$-$\mathrm{N_2}$')

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

# fig.savefig('output/plots/multiGasEuler/convEuler.pdf')
