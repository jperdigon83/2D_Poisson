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

test_data = ascii.read('result.dat',format='no_header')
m= np.int(np.sqrt(len(test_data)))

xx =np.reshape(test_data['col1'],(m,m))
yy =np.reshape(test_data['col2'],(m,m))
u = np.reshape(test_data['col3'],(m,m))

Lx = xx[-1,0] + (xx[-1,0]-xx[-2,0])
Ly = yy[0,-1] + (yy[0,-1]-yy[0,-2])
uth = np.zeros((m,m))
K = 20

for i in range(m):
    for j in range(m):

        for k in range(K):
            uth[i,j] += np.sin((2*k+1)*np.pi*0.5*(2*Lx*xx[i,j]))*(np.sinh((2*k+1)*np.pi*0.5*(2*Ly*yy[i,j])) + np.sinh((2*k+1)*np.pi*0.5*(2-2*Ly*yy[i,j])) )/((2*k+1)**3*np.sinh(np.pi*(2*k+1)))
        
        uth[i,j] = 0.5*(1-(2*xx[i,j]-1)**2) - 16*uth[i,j]/np.pi**3


plt.figure(figsize=(12,5))

ax1 = plt.subplot(121)
ax1.pcolor(xx,yy,u)

ax2 = plt.subplot(122)
pcolor(xx,yy,uth)

plt.tight_layout()





