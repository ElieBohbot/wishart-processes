#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 18:27:46 2017

@author: EB
"""

import pylab as py
from WishartDim2 import *
import numpy as np
import matplotlib.pyplot as plt


T = 1.0
N = 1
a = 4.5


#matrice v fixee    
x = np.array([[0.04,0.02],[0.02,0.04]])
u = 0.2
v = 0.5  

sig = WIS2_array_I2(a,x,T,N)

print(sig)    