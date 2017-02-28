#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------ WishartDim2_O2_e1_Test2.py -------------------

Created on Fri Feb 24 23:40:47 2017

@author: EB
"""

import pylab as py
from WishartDim2_O2 import *
from WishartDim2 import *
import numpy as np
import matplotlib.pyplot as plt



M = 500

#Parametres de test
T = 1.0
a = 2.5
x = np.array([[2,1],[1,1]])

#matrice z fixee, parametres u et v    
u = 0.2
v = 0.5 

#Fonction f
#A FAIRE : Calcul exact de E[f(X_T)] pour une fonction f donnee et ces parametres
def MaFonction(X):
    #return py.det(X)
    return np.cos(X[0][0])

moyenneEmp = np.zeros(M)

for i in np.arange(M):
    X = WIS2_StepByStep(a,x,0,T,u,v)
    moyenneEmp[i] = MaFonction(X)
    
mean = np.mean(moyenneEmp)
integers1toM = np.arange(1,M+1)
moyenneEmp = np.cumsum(moyenneEmp)/integers1toM

#nombre de points de discretisation de [0,1]
N_MAX = 2000
pas = 10
N = N_MAX/pas
moyenneEmp_O2 = np.zeros(N)
integers1toN = np.arange(1,N+1)

j=0
for n in np.arange(1,N_MAX,pas):
    sig = WIS2_O2_array(a,x,T,n,u,v)
    moyenneEmp_O2[j] = MaFonction(sig[n])
    j+=1

mean_O2 = np.mean(moyenneEmp_O2)
moyenneEmp_O2 = np.cumsum(moyenneEmp_O2)/integers1toN

# Affichage des trajectoires
#plt.subplot(211)
#plt.plot(integers1toM, moyenneEmp, lw=2, color="b")
#plt.axhline(mean, color="r", label="det_emp_simul_exact")
#plt.legend(loc="best")
#
#plt.subplot(212)
plt.plot(integers1toN, moyenneEmp_O2, color="g")
plt.axhline(mean, color="r", label="moyEmp_simul_exact")
plt.legend(loc="best")

plt.show()








