#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 20:03:27 2017

@author: EB
"""

import numpy as np

#Dans ce cas k = 0

sig = 2.0

def Y():
    U = np.random.rand()
    return (U < 1./6)*(-np.sqrt(3)) + ((U>1./6)&(U<1./3))*np.sqrt(3)
    
def phi(a,x,t,w):
    sqr = np.sqrt((a-sig**2/4)*t/2 + x)
    sqr = sqr + sig/2*w
    return sqr**2 + (a-sig**2/4)*t/2

def K2(a,t):
    if (sig**2>4*a):
        sqr = np.sqrt((sig**2/4 - a)*t/2)
        sqr = sqr + sig/2*np.sqrt(3*t)
        sqr = sqr**2
        return (sig**2/4 - a)*t/2 + sqr
    else:
        return 0
   
def u1(a,t,x):
    return x + a*t
    
def u2(a,t,x):
    return u1(a,t,x)**2 + sig**2*t*(a*t/2+x)
    
def pi(a,t,x):
    delta = 1 - u1(a,t,x)**2/u2(a,t,x)
    return (1 - np.sqrt(delta))/2

def CIR_O2_StepByStep(a,t,x):
    if (x >= K2(a,t)):
        Z = Y()
        return phi(a,x,t,np.sqrt(t)*Z)
    else:
        U = np.random.rand()
        return (U < pi(a,t,x))*u1(a,t,x)/(2*pi(a,t,x)) + (U > pi(a,t,x))*u1(a,t,x)/(2*(1 - pi(a,t,x)))
    
    
