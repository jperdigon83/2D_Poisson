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


plt.close('all')

serial_error = ascii.read("serial_error.dat",format='commented_header')
serial_time = ascii.read("time_serial.dat",format='commented_header')

N = serial_error['N']
L2err = serial_error['errL2']
t = serial_time['t']

plt.figure(figsize=(17,8))

ax1=plt.subplot(121)

plt.scatter(N,L2err)
plt.plot(N,L2err[0]*N[0]**2/N**2,'--k',label=r'$\propto n^{-2}$')
plt.xscale('log',basex=2)
plt.yscale('log')
plt.xlabel(r'$n$',fontsize=16)
plt.ylabel(r'$L_2 \ error$',fontsize=16)
plt.legend(fontsize=16)
ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))


ax2=plt.subplot(122)

plt.scatter(N,t)
plt.plot(N,t[-1]*(N/N[-1])**2*np.log2(N)/np.log2(N[-1]),'--k',label=r'$\propto n^2 \log(n)$')
plt.xscale('log',basex=2)
plt.yscale('log')
plt.xlabel(r'$n$',fontsize=16)
plt.ylabel(r'$CPU-time \ ( s)$',fontsize=16)
plt.tight_layout()
plt.legend(fontsize=16)
ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))


