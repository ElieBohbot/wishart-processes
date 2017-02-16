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
N = 10

a = 5
u = 1
v = 2

r = 0.03
sigma = 0.4
rho = 0
x = np.array([[sigma**2,-rho],[-rho,sigma**2]])
S_0 = 100

sig = WIS2_array(a,x,T,N,u,v)

print(lin.sqrtm(T/N*sig[9]))
#%%

S = np.zeros((2,N))
S[0][0] = S_0
S[1][0] = S_0

for i in 1+np.arange(N-1):
    X = lin.sqrtm(T/N*sig[i]) #sig = XX'
    G = np.random.randn(2)
    G = np.dot(X,G)
    S[:,i] = S[:,i-1]*np.exp(r*T/N*np.ones(2) + G - 0.5*T/N*np.diag(sig[i]))
    
integers1toN = np.arange(1,N+1)
t = np.linspace(0,T,N) 

plt.plot(t,S[0])
plt.plot(t,S[1])
plt.grid()

