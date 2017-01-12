#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 18:27:46 2017

@author: EB
"""

import pylab as py
from WishartDim2 import *
import numpy as np

M = 5000

reals = np.zeros(M)

def expoTrace(v,X):
    return np.exp(py.trace(py.dot(v,X)))

v = np.array([[0.3,0.2],[0.2,0.4]])   
    
for i in np.arange(M):
    res = WIS2_array_e1(a,x,T,N)
    reals[i] = expoTrace(res[N],v)
    
val_approx = np.mean(reals)
val_exacte = LaplaceTransform(v,T,a,x)

print(val_approx)
print(val_exacte)
    