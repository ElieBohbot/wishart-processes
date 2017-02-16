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


T = 0.01
N = 1
a = 2.5
x = np.array([[2,1],[1,1]])

M = 200
N = 50

val_approx = np.zeros(N)
std_approx = np.zeros(N)

#Fonction caracteristique cas a = u*I2 et b = v*I2
def LaplaceTransform(z,t,a,x,u,v):
    b = v*np.eye(2)
    q = np.array([[u**2/(2*v)*(np.exp(2*v*t)-1),0],
                  [0,u**2/(2*v)*(np.exp(2*v*t)-1)]])
    m = np.exp(b*t)
    e = np.eye(2)-2*py.dot(q,z)
    w = py.inv(e)
    w = np.dot(w,np.dot(m,np.dot(x,m)))
    y = py.dot(z,w)
    den = py.det(e)**(a/2.0)
    num = np.exp(py.trace(y))
    return num/den


#Fonction caracteristique empirique
def expoTrace(v,X):
    return np.exp(py.trace(py.dot(v,X)))

#matrice v fixee    
z = np.array([[0.03,0.02],[0.02,0.04]])
u = 0.2
v = 0.5   
    
def SimulTransform(lam,i):
    reals = np.zeros(M)
    for j in np.arange(M):
       res = WIS2_StepByStep(a,x,0,T,u,v)
       reals[j] = expoTrace(res,lam[i]*z)
    val_approx[i]=np.mean(reals)
    std_approx[i]=np.std(reals)


lam = np.linspace(2,3,N)

for i in np.arange(N):
    SimulTransform(lam,i)

def LaplaceTransformFunc(lam):
    res = np.zeros(len(lam))
    for i in np.arange(len(lam)):
        res[i]=LaplaceTransform(lam[i]*z,T,a,x,u,v)
    return res
  
#demi largeur et intervalle de confiance    
demi_largeur = 1.96*std_approx/np.sqrt(M)
IC_sup = val_approx + demi_largeur
IC_inf = val_approx - demi_largeur


#Affichage    
plt.plot(lam, LaplaceTransformFunc(lam), color='blue', label='Laplace Transform')
plt.plot(lam, val_approx, color='green', label='Laplace Transform emp')
plt.plot(lam, IC_sup, color='red', label='IC a 95%')
plt.plot(lam, IC_inf, color='red')
plt.legend(loc='best')
plt.show()

    