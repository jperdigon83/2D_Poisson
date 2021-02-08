# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:49:51 2021

@author: jperdigon
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from pylab import *
import matplotlib.ticker as mtick
from mpl_toolkits.mplot3d import Axes3D


plt.close('all')


# =============================================================================
# 
# =============================================================================
"""
plt.figure(figsize=(17,8))

N = ascii.read('time_parallel_P=1.dat',format='commented_header')['N']

t1 = ascii.read('time_parallel_P=1.dat',format='commented_header')['t']
t2 = ascii.read('time_parallel_P=2.dat',format='commented_header')['t']
t4 = ascii.read('time_parallel_P=4.dat',format='commented_header')['t']
t6 = ascii.read('time_parallel_P=6.dat',format='commented_header')['t']

P = np.array([1,2,4,6],dtype='int')
Pid = np.linspace(1,8,128)

TP_N32 = np.array([1,t1[-5]/t2[-5],t1[-5]/t4[-5],t1[-5]/t6[-5]])
TP_N64 = np.array([1,t1[-4]/t2[-4],t1[-4]/t4[-4],t1[-4]/t6[-4]])
TP_N128 = np.array([1,t1[-3]/t2[-3],t1[-3]/t4[-3],t1[-3]/t6[-3]])
TP_N256 = np.array([1,t1[-2]/t2[-2],t1[-2]/t4[-2],t1[-2]/t6[-2]])
TP_N512 = np.array([1,t1[-1]/t2[-1],t1[-1]/t4[-1],t1[-1]/t6[-1]])

ax1 = plt.subplot(111)

plt.scatter(N,t1,label=r'$P=1$')
plt.plot(N,t1[-1]*(N/N[-1])**2*np.log2(N)/np.log2(N[-1]),'--k',label=r'$\propto n^2 \log(n)$')
plt.scatter(N,t2,label=r'$P=2$')
plt.plot(N,t2[-1]*(N/N[-1])**2*np.log2(N)/np.log2(N[-1]),'--k')
plt.scatter(N,t4,label=r'$P=4$')
plt.plot(N,t4[-1]*(N/N[-1])**2*np.log2(N)/np.log2(N[-1]),'--k')
plt.scatter(N,t6,label=r'$P=6$')
plt.plot(N,t6[-1]*(N/N[-1])**2*np.log2(N)/np.log2(N[-1]),'--k')

plt.xscale('log',basex=2)
plt.xlabel(r'$n$',fontsize=16)
ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))

plt.yscale('log')
plt.ylabel(r'$T_P \ (s)$',fontsize=16)

plt.legend(fontsize=16)


plt.figure(figsize=(17,8))

ax1 = plt.subplot(121)

ax1.scatter(P,TP_N32,label=r'$n=32$')
ax1.scatter(P,TP_N64,label=r'$n=64$')
ax1.scatter(P,TP_N128,label=r'$n=128$')
ax1.scatter(P,TP_N256,label=r'$n=256$')
plt.scatter(P,TP_N512,label=r'$n=512$')
ax1.plot(Pid,Pid,'k--',label=r'$ideal \ , \ T_P = P$')

plt.xscale('log',basex=2)
plt.xlabel(r'$P$',fontsize=16)
ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))

plt.ylabel(r'$ T_1 \ / \ T_P$',fontsize=16)

plt.legend(fontsize=16)

ax2 = plt.subplot(122,sharey=ax1)

P = np.linspace(1,32,200)

for i in range(2,len(N)):
    
    tid = N[i]**2*np.log(N[i])
    tr = N[i]**2*np.log(N[i])/P + 1e-1*(P-1)*N[i]**2
    
    plt.plot(P,tid/tr)
plt.plot(P,P,'k--')
plt.xscale('log',basex=2)
#plt.yscale('log',basey=2)
ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
plt.xlabel(r'$P$',fontsize=16)
plt.tight_layout()

"""

# =============================================================================
# # Plot of the error and computing time wrt n, for P=1

# =============================================================================
"""
L2err =  ascii.read('error_parallel_P=1.dat',format='commented_header')['errL2']

plt.figure(figsize=(17,8))

ax1=plt.subplot(121)

plt.scatter(N,L2err)
plt.plot(N,L2err[0]*N[0]**2/N**2,'--k',label=r'$\propto n^{-2}$')
plt.xscale('log',basex=2)
plt.yscale('log')
plt.xlabel(r'$n$',fontsize=16)
plt.ylabel(r'$\epsilon_2$',fontsize=16)
plt.legend(fontsize=16)
ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))

ax2=plt.subplot(122)

plt.scatter(N,t1)
plt.plot(N,t1[-1]*(N/N[-1])**2*np.log2(N)/np.log2(N[-1]),'--k',label=r'$\propto n^2 \log(n)$')
plt.xscale('log',basex=2)
plt.yscale('log')
plt.xlabel(r'$n$',fontsize=16)
plt.ylabel(r'$T_1 \ ( s)$',fontsize=16)
plt.tight_layout()
plt.legend(fontsize=16)
ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
"""

# =============================================================================
# 
# =============================================================================
dataU1 = ascii.read('res_f=1.dat',format="commented_header")
N = np.int(np.sqrt(len(dataU1)))
xx = np.reshape(dataU1['x'],(N,N))
yy = np.reshape(dataU1['y'],(N,N))
u1 = np.reshape(dataU1['u'],(N,N))

dataU2 = ascii.read('res_f=sinsin.dat',format="commented_header")
u2 = np.reshape(dataU2['u'],(N,N))

dataU3 = ascii.read('res_f=expsinsin.dat',format="commented_header")
u3 = np.reshape(dataU3['u'],(N,N))

dataU4 = ascii.read('res_f=pointsource.dat',format="commented_header")
u4 = np.reshape(dataU4['u'],(N,N))


fig = plt.figure(figsize=(17,8))

ax1 = fig.add_subplot(221, projection='3d')
ax1.plot_surface(xx, yy, u1/np.max(u1), cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax1.set_xlabel(r'$x$',fontsize=16)
ax1.set_ylabel(r'$y$',fontsize=16)
ax1.set_zlabel(r'$\frac{u}{max(u)}$',fontsize=16)

ax2 = fig.add_subplot(222, projection='3d')
ax2.plot_surface(xx, yy, u2/np.max(u2), cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax2.set_xlabel(r'$x$',fontsize=16)
ax2.set_ylabel(r'$y$',fontsize=16)
ax2.set_zlabel(r'$\frac{u}{max(u)}$',fontsize=16)


ax3 = fig.add_subplot(223, projection='3d')
ax3.plot_surface(xx, yy, u3/np.max(u3), cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax3.set_xlabel(r'$x$',fontsize=16)
ax3.set_ylabel(r'$y$',fontsize=16)
ax3.set_zlabel(r'$\frac{u}{max(u)}$',fontsize=16)


ax4 = fig.add_subplot(224, projection='3d')
ax4.plot_surface(xx, yy, u4/np.max(u4), cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax4.set_xlabel(r'$x$',fontsize=16)
ax4.set_ylabel(r'$y$',fontsize=16)
ax4.set_zlabel(r'$\frac{u}{max(u)}$',fontsize=16)




plt.tight_layout()


