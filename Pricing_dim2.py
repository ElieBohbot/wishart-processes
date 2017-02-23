#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 20:27:43 2017

@author: EB
"""

import pylab as py
import numpy as np
import scipy.linalg as lin
from WishartDim2 import *
import matplotlib.pyplot as plt

T = 1.0
N = 30

a = 4.5
x = np.array([[0.04,0.02],[0.02,0.04]])
u = 0.2
v = 0.5

r = 0.03
sigma = np.sqrt(0.1)
rho = 0

S_0 = 100

sig = WIS2_array(a,x,T,N,u,v)
#sig = WIS2_array_e1(a,x,T,N)
#sig = WIS2_array_I2(a,x,T,N)


S = np.zeros((2,N))
S[0][0] = S_0
S[1][0] = S_0

for i in 1+np.arange(N-1):
    X = lin.sqrtm(T/N*sig[i]) #sig = XX'
    G = np.random.randn(2)
    G = np.dot(X,G)
    S[:,i] = S[:,i-1]*np.exp(r*T/N*np.ones(2) + G - 0.5*T/N*np.diag(sig[i]))
    
t = np.linspace(0,T,N) 

plt.plot(t,S[0])
plt.plot(t,S[1])
plt.grid()

