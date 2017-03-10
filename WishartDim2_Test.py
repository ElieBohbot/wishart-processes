#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------ WishartDim2_Test.py -------------------

Created on Thu Jan 12 18:27:46 2017

@author: EB
"""

import pylab as py
from WishartDim2 import *
import numpy as np
import matplotlib.pyplot as plt
from time import time

#==============================================================================
# Test de la simulation WIS_2(x,a,b,A;t) ou A=u*I_2 et b=v*I_2
# --------------------------------------------------------------
# On discretise l'intervalle [2,3] avec N points 
# Pour chaque point lambda de cet intervalle, on calcule la fonction caracteristique
# au point lambda*z ou z est une matrice fixee avec des coefficients assez faibles 
# pour eviter une explosion des valeurs de la fonction caracteristique (a cause de l'exponentielle)
# On calcule egalement la transformee de Laplace empirique en ce point par une methode de 
# Monte-Carlo avec M simulations.
# On affiche la transformee de Laplace et sa fonction empirique ainsi que les intervalles 
# de confiance sur [2,3]    
#==============================================================================

#Parametres de test
T = 1.0
a = 2.5
x = np.array([[1,0],[0,1]])
M = int(1.e4)
N = 20

#matrice z fixee, parametres u et v    
z = np.array([[0.03,0.02],[0.02,0.04]])
u = 0.2
v = 0.5 

val_approx = np.zeros(N)
std_approx = np.zeros(N)
reals = np.zeros(M)

#Fonction caracteristique cas a = u*I2 et b = v*I2
def LaplaceTransform(z,t,a,x,u,v):
    q = np.array([[u**2/(2*v)*(np.exp(2*v*t)-1),0],
                  [0,u**2/(2*v)*(np.exp(2*v*t)-1)]])
    m = np.array([[np.exp(t*v),0],
                  [0,np.exp(t*v)]])
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

 
#Simulation des M Wishart
simul = np.zeros((M,2,2)) 
for i in np.arange(M):
    simul[i]=WIS2_StepByStep(a,x,0,T,u,v)

    
def SimulTransform(lam,i):
    for j in np.arange(M):
       reals[j] = expoTrace(lam[i]*z,simul[j])
    val_approx[i]= np.mean(reals)
    std_approx[i]= np.std(reals)


lam = np.linspace(6,8,N)

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

    