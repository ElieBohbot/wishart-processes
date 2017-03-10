#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------ WishartDim2.py -------------------

Created on Mon Jan  2 18:25:19 2017

@author: EB
"""

import pylab as py
import numpy as np
from SQB import *

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
        U[0][0]=X[0][0]-(X[0][1]**2)/X[1][1]
        U[0][1]=X[0][1]/np.sqrt(X[1][1])
        U[1][1]=X[1][1]
        #
        temp = np.zeros((2,2))
        temp[0][0] = BESQ(a-1,U[0][0],t_1,t_2)
        temp[0][1] = U[0][1] + np.sqrt(t_2-t_1)*np.random.randn()
        temp[1][1]= U[1][1]
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

#==============================================================================
# Algorithm 2 : Exact simulation WIS_2(x,a,0,I_2;t)
#==============================================================================  
    
def WIS2_StepByStep_I2(a,X,t_1,t_2):
    y = X
    y = WIS2_StepByStep_e1(a,y,t_1,t_2)
    p = np.array([[0,1],[1,0]])
    Y = WIS2_StepByStep_e1(a,py.dot(py.dot(p,y),p),t_1,t_2)
    y = py.dot(py.dot(p,Y),p)
    return y
    
def WIS2_array_I2(a,x,T,N):
    t = np.linspace(0,T,N+1)
    res = np.zeros((N+1,2,2))
    res[0] = x
    for i in np.arange(N):
        res[i+1] = WIS2_StepByStep_I2(a,res[i],t[i],t[i+1])
    return res
    
    
#==============================================================================
# Algorithm 3 : Exact simulation WIS_2(x,a,b,A;t)
#============================================================================== 

    
#on suppose que A=u*I_2 et b=v*I_2  
def WIS2_StepByStep(a,X,t_1,t_2,u,v):
    #Calcul de q_t
    t = t_2 - t_1
    q = np.array([[u**2/(2*v)*(np.exp(2*v*t)-1),0],
                  [0,u**2/(2*v)*(np.exp(2*v*t)-1)]])
    c = np.sqrt(q/t)
    m = np.array([[np.exp(t*v),0],
                  [0,np.exp(t*v)]])
    x = np.dot(m,np.dot(X,m))
    c_inv = py.inv(c)
    x_prime = np.dot(c_inv,np.dot(x,c_inv))
    Y = WIS2_StepByStep_I2(a,x_prime,t_1,t_2)
    return np.dot(c,np.dot(Y,c))
    

def WIS2_array(a,x,T,N,u,v):
    t = np.linspace(0,T,N+1)
    res = np.zeros((N+1,2,2))
    res[0] = x
    for i in np.arange(N):
        res[i+1] = WIS2_StepByStep(a,res[i],t[i],t[i+1],u,v)
    return res       


    

    
    


    
    
    
    
    
    

