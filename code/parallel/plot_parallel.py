# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:49:51 2021

@author: jperdigon
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from pylab import *

plt.close('all')

# test_data = ascii.read('res.dat',format='commented_header')
# m= np.int(np.sqrt(len(test_data)))

# xx =np.reshape(test_data['col1'],(m,m))
# yy =np.reshape(test_data['col2'],(m,m))
# u = np.reshape(test_data['col3'],(m,m))

# dx = xx[1,0] - xx[0,0]
# dy = yy[0,1] - yy[0,0]

# M = 20
# N = 20
# uth = np.zeros((m,m))

# for o in range(1,M):
#     for n in range(1,N):
#         coeff = -4*(1-(-1)**o)*(1-(-1)**n)/((o**2+n**2)*o*n*np.pi**2)
#         uth += coeff * np.sin(o*xx) * np.sin(n*yy)

# diff = np.zeros((m,m))
# diff = np.abs( (uth - u)/uth)

# plt.figure(figsize=(12,5))

# ax1 = plt.subplot(121)
# ax1.pcolor(xx,yy,u)

# ax2 = plt.subplot(122)
# pcolor(xx,yy,uth)
# plt.tight_layout()

# print( np.sqrt(np.sum((u-uth)**2)/(m*m) )  )

err_data = ascii.read("errL2.dat",format='commented_header')
N = err_data['N']
L2err = err_data['errL2']
t = err_data['t']

plt.figure(figsize=(17,8))

plt.subplot(121)

plt.scatter(N,L2err)
plt.plot(N,L2err[0]*N[0]**2/N**2,'--k')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$N$',fontsize=16)
plt.ylabel(r'$L_2 \ error$',fontsize=16)

plt.subplot(122)

plt.scatter(N,t/1000)
plt.plot(N,t[-1]*(N/N[-1])**2*np.log2(N)/np.log2(N[-1])*1./1000,'--k')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$N$',fontsize=16)
plt.ylabel(r'$CPU-time \ ( s)$',fontsize=16)
plt.tight_layout()


