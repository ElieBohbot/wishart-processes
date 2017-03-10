#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------ WishartDim2_O2_I2_Test.py -------------------

Created on Thu Jan 12 18:27:46 2017

@author: EB
"""

import pylab as py
from WishartDim2 import *
import numpy as np
import matplotlib.pyplot as plt
from time import time

#==============================================================================
# Test de la simulation WIS_2(x,a,0,I_2;t)
# ------------------------------------------
# On discretise l'intervalle [2,3] avec N points 
# Pour chaque point lambda de cet intervalle, on calcule la fonction caracteristique
# au point lambda*v ou v est une matrice fixee avec des coefficients assez faibles 
# pour eviter une explosion des valeurs de la fonction caracteristique (a cause de l'exponentielle)
# On calcule egalement la transformee de Laplace empirique en ce point par une methode de 
# Monte-Carlo avec M simulations.
# On affiche la transformee de Laplace et sa fonction empirique ainsi que les intervalles 
# de confiance sur [2,3]    
#==============================================================================

t_0 = time()

#Discretisation de [2,3] (N) et nombre de simulations par point (M)
M = int(1.e4)
N = 20
N_array = 10

#Parametres de test
T = 1.0
a = 1.1
x = np.array([[1,0],[0,1]])

#matrice v fixee
v = np.array([[0.04,0.02],[0.02,0.04]])  

val_approx = np.zeros(N)
std_approx = np.zeros(N)
reals = np.zeros(M)

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
v = np.array([[0.03,0.02],[0.02,0.04]])   
    
#Simulation des M Wishart
simul = np.zeros((M,N_array,2,2)) 
for i in np.arange(M):
    simul[i]=WIS2_O2_array_I2(a,x,T,N_array-1)
    

def SimulTransform(lam,i):
    for j in np.arange(M):
       reals[j] = expoTrace(lam[i]*v,simul[j][N_array-1])
    val_approx[i]= np.mean(reals)
    std_approx[i]= np.std(reals)


lam = np.linspace(0,5,N)

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
plt.plot(lam, val_approx, color='green', label='Laplace Transform emp_O2 I2')
plt.plot(lam, IC_sup, color='red', label='IC a 95%')
plt.plot(lam, IC_inf, color='red')
plt.legend(loc='best')
plt.show()

t_1 = time()
t = t_1 - t_0
print(t)