# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:49:51 2021

@author: jperdigon
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii



m = 128

test_data = ascii.read('sol_f=1.dat',format='no_header')['col1']
u = np.reshape(test_data,(m,m))

plt.imshow(u)