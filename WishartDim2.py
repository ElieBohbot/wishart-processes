#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 18:25:19 2017

@author: EB
"""

import pylab as py
import numpy as np
from SQB import *

T = 0.01
N = 50
a = 1.5
x = np.array([[2,1],[1,1]])

#==============================================================================
# Algorithm 1 : Exact simulation WIS_2(x,a,0,e^1_2;t)
#==============================================================================

def WIS2_StepByStep_e1(a,X,t_1,t_2):
    res = np.zeros((2,2))
    if X[1][1]==0:
        res[1][1]=X[1][1]
        res[1][0]=X[1][0]
        res[0][1]=X[0][1]
        res[0][0]=BESQ(a,X[0][0],t_1,t_2)
    else:
        U = np.zeros((2,2))
        U[0][0]=X[0][0]-X[0][1]**2/X[1][1]
        U[0][1]=X[0][1]/np.sqrt(X[1][1])
        U[1][1]=X[1][1]
        #
        temp = np.zeros((2,2))
        temp[0][0] = BESQ(a-1,U[0][0],t_1,t_2)
        temp[0][1] = U[0][1] + np.sqrt(t_2-t_1)*np.random.randn()
        temp[1][1]=U[1][1]
        #
        res[0][0] = temp[0][0]+temp[0][1]**2
        res[1][1] = temp[1][1]
        res[0][1] = temp[0][1]*np.sqrt(temp[1][1])
        res[1][0] = temp[0][1]*np.sqrt(temp[1][1])
    return res
    

def WIS2_array_e1(a,x,T,N):
    t = np.linspace(0,T,N+1)
    res = np.zeros((N+1,2,2))
    res[0] = x
    for i in np.arange(N):
        res[i+1] = WIS2_StepByStep_e1(a,res[i],t[i],t[i+1])
    return res

def WIS2_StepByStep_I2(a,X,t_1,t_2):
    y = X
    y = WIS2_StepByStep_e1(a,y,t_1,t_2)
    p = np.array([[0,1],[1,0]])
    Y = WIS2_StepByStep_e1(a,py.dot(py.dot(p,y),p),t_1,t_2)
    y = py.dot(py.dot(p,Y),p)
    return y
    
def WIS2_StepByStep(a,X,t_1,t_2,A,b):
    
        

def LaplaceTransform(v,t,a,x):
    q = np.array([[t,0],[0,0]])
    z = np.eye(2)-2*py.dot(q,v)
    w = py.inv(z)
    y = py.dot(v,w)
    den = py.det(z)**(a/2.0)
    num = np.exp(py.trace(py.dot(y,x)))
    return num/den
    

    
    


    
    
    
    
    
    

