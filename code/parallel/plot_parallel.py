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

N = ascii.read('time_parallel_P=1.dat',format='commented_header')['N']

t1 = ascii.read('time_parallel_P=1.dat',format='commented_header')['t']
t2 = ascii.read('time_parallel_P=2.dat',format='commented_header')['t']
t4 = ascii.read('time_parallel_P=4.dat',format='commented_header')['t']
t6 = ascii.read('time_parallel_P=6.dat',format='commented_header')['t']

P = np.array([1,2,4,6],dtype='int')
Pid = np.linspace(1,6,128)

TP_N32 = np.array([1,t1[-5]/t2[-5],t1[-5]/t4[-5],t1[-5]/t6[-5]])
TP_N64 = np.array([1,t1[-4]/t2[-4],t1[-4]/t4[-4],t1[-4]/t6[-4]])
TP_N128 = np.array([1,t1[-3]/t2[-3],t1[-3]/t4[-3],t1[-3]/t6[-3]])
TP_N256 = np.array([1,t1[-2]/t2[-2],t1[-2]/t4[-2],t1[-2]/t6[-2]])
TP_N512 = np.array([1,t1[-1]/t2[-1],t1[-1]/t4[-1],t1[-1]/t6[-1]])

plt.figure(figsize=(17,8))

ax1 = plt.subplot(121)

plt.plot(N,t1,label='$P=1$')
plt.plot(N,t2,label='$P=2$')
plt.plot(N,t4,label='$P=4$')
plt.plot(N,t6,label='$P=6$')

plt.xscale('log',basex=2)
plt.xlabel(r'$n$',fontsize=16)
ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))

plt.yscale('log')
plt.ylabel(r'$T_P \ (s)$',fontsize=16)

plt.legend(fontsize=16)

ax2 = plt.subplot(122)

plt.plot(P,TP_N32,label='n=32')
plt.plot(P,TP_N64,label='n=64')
plt.plot(P,TP_N128,label='n=128')
plt.plot(P,TP_N256,label='n=256')
plt.plot(P,TP_N512,label='n=512')
plt.plot(Pid,Pid,'k--')

plt.xscale('log',basex=2)
plt.xlabel(r'$P$',fontsize=16)
ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))

plt.ylabel(r'$ T_1 \ / \ T_P$',fontsize=16)

plt.legend(fontsize=16)
plt.tight_layout()
