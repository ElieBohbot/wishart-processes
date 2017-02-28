#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------ WishartDim2_O2_e1_Test.py -------------------

Created on Thu Jan 12 18:27:46 2017

@author: EB
"""

import pylab as py
from WishartDim2_O2 import *
import numpy as np
import matplotlib.pyplot as plt
from time import time

#==============================================================================
# Test de la simulation WIS_2(x,a,0,e^1_2;t) avec un schema a l'ordre 2
# ---------------------------------------------------------------------
# On discretise l'intervalle [2,3] avec N points 
# Pour chaque point lambda de cet intervalle, on calcule la fonction caracteristique
# au point lambda*v ou v est une matrice fixee avec des coefficients assez faibles 
# pour eviter une explosion des valeurs de la fonction caracteristique (a cause de l'exponentielle)
# On calcule egalement la transformee de Laplace empirique en ce point par une methode de 
# Monte-Carlo avec M simulations.
# On affiche la transformee de Laplace et sa fonction empirique ainsi que les intervalles 
# de confiance sur [2,3]    
#
# Contrairement a la simulation exacte on constate qu'on n'observe pas de convergence en loi du 
# Wishart simule, ce qui est en accord avec les resultats theoriques
#==============================================================================

t_0 = time()

#Discretisation de [2,3] (N) et nombre de simulations par point (M)
M = 5000
N = 30

T = 1.0
a = 2.5
x = np.array([[2,1],[1,1]])

val_approx = np.zeros(N)
std_approx = np.zeros(N)

#Fonction caracteristique cas a = e1
def LaplaceTransform_e1(v,t,a,x):
    q = np.array([[t,0],[0,0]])
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
    
def SimulTransform(lam,i):
    reals = np.zeros(M)
    for j in np.arange(M):
       res = WIS2_O2_StepByStep_e1(a,x,0,T)
       reals[j] = expoTrace(res,lam[i]*v)
    val_approx[i]=np.mean(reals)
    std_approx[i]=np.std(reals)


lam = np.linspace(2,3,N)

for i in np.arange(N):
    SimulTransform(lam,i)

def LaplaceTransformFunc(lam):
    res = np.zeros(len(lam))
    for i in np.arange(len(lam)):
        res[i]=LaplaceTransform_e1(lam[i]*v,T,a,x)
    return res
  
#demi largeur et intervalle de confiance    
demi_largeur = 1.96*std_approx/np.sqrt(M)
IC_sup = val_approx + demi_largeur
IC_inf = val_approx - demi_largeur


#Affichage    
plt.plot(lam, LaplaceTransformFunc(lam), color='blue', label='Laplace Transform e1')
plt.plot(lam, val_approx, color='green', label='Laplace Transform emp_O2 e1')
plt.plot(lam, IC_sup, color='red', label='IC a 95%')
plt.plot(lam, IC_inf, color='red')
plt.legend(loc='best')
plt.show()

t_1 = time()
t = t_1 - t_0
print(t)
