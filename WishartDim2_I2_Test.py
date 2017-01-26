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

M = 200
N = 50

val_approx = np.zeros(N)
std_approx = np.zeros(N)

#Fonction caracteristique cas a=I2 et b=0
def LaplaceTransform_I2(v,t,a,x):
    q = np.array([[t,0],[0,t]])
    z = np.eye(2)-2*py.dot(q,v)
    w = py.inv(z)
    y = py.dot(v,w)
    den = py.det(z)**(a/2.0)
    num = np.exp(py.trace(py.dot(y,x)))
    return num/den


#Fonction caracteristique empirique
def expoTrace(v,X):
    return np.exp(py.trace(py.dot(v,X)))

#matrice v fixee    
v = np.array([[0.3,0.2],[0.2,0.4]])   
    
def SimulTransform(lam,i):
    reals = np.zeros(M)
    for j in np.arange(M):
       res = WIS2_StepByStep_I2(a,x,0,T)
       reals[j] = expoTrace(res,lam[i]*v)
    val_approx[i]=np.mean(reals)
    std_approx[i]=np.std(reals)


lam = np.linspace(2,3,N)

for i in np.arange(N):
    SimulTransform(lam,i)

def LaplaceTransformFunc(lam):
    res = np.zeros(len(lam))
    for i in np.arange(len(lam)):
        res[i]=LaplaceTransform_I2(lam[i]*v,T,a,x)
    return res
  
#demi largeur et intervalle de confiance    
demi_largeur = 1.96*std_approx/np.sqrt(M)
IC_sup = val_approx + demi_largeur
IC_inf = val_approx - demi_largeur


#Affichage    
plt.plot(lam, LaplaceTransformFunc(lam), color='blue', label='Laplace Transform I2')
plt.plot(lam, val_approx, color='green', label='Laplace Transform emp I2')
plt.plot(lam, IC_sup, color='red', label='IC a 95%')
plt.plot(lam, IC_inf, color='red')
plt.legend(loc='best')
plt.show()

    