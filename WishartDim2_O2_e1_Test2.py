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


#Parametres de test
T = 1.0
a = 2.5
x = np.array([[2,1],[1,1]])

#Fonction f
#A FAIRE : Calcul exact de E[f(X_T)] pour une fonction f donnee et ces parametres
def MaFonction(X):
    return py.det(X)
    #return py.trace(X)
    #return np.cos(X[0][0])

M = 5000
moyenneEmp = np.zeros(M)

for i in np.arange(M):
    X = WIS2_StepByStep_e1(a,x,0,T)
    moyenneEmp[i] = MaFonction(X)
    
mean = np.mean(moyenneEmp)
integers1toM = np.arange(1,M+1)
moyenneEmp = np.cumsum(moyenneEmp)/integers1toM

#demi largeur et intervalle de confiance    
std_approx = np.std(moyenneEmp)
demi_largeur = 1.96*std_approx/np.sqrt(M)
IC_sup = moyenneEmp + demi_largeur
IC_inf = moyenneEmp - demi_largeur

#nombre de points de discretisation de [0,1]
N_MAX = 200

M_O2 = 2000
moyenneEmp_O2 = np.zeros(M_O2)
integers1toM_O2 = np.arange(1,M_O2+1)

for i in np.arange(M_O2):
    sig = WIS2_O2_array_e1(a,x,T,N_MAX)
    moyenneEmp_O2[i] = MaFonction(sig[N_MAX])


mean_O2 = np.mean(moyenneEmp_O2)
moyenneEmp_O2 = np.cumsum(moyenneEmp_O2)/integers1toM_O2

#demi largeur et intervalle de confiance    
std_approx_O2 = np.std(moyenneEmp_O2)
demi_largeur_O2 = 1.96*std_approx_O2/np.sqrt(M_O2)
IC_sup_O2 = moyenneEmp_O2 + demi_largeur_O2
IC_inf_O2 = moyenneEmp_O2 - demi_largeur_O2

# Affichage des trajectoires
plt.figure()
plt.plot(integers1toM, moyenneEmp, lw=2, color="b")
#plt.plot(integers1toM, IC_sup, color='red', label='IC a 95%')
#plt.plot(integers1toM, IC_inf, color='red')
#plt.axhline(mean, color="r", label="det_emp_simul_exact")
plt.legend(loc="best")
#lt.show()
#
plt.figure()
plt.plot(integers1toM_O2, moyenneEmp_O2, lw=1, color="g")
#plt.plot(integers1toM_O2, IC_sup_O2, color='red', label='IC a 95%')
#plt.plot(integers1toM_O2, IC_inf_O2, color='red')
plt.axhline(mean, color="r", label="moyEmp_simul_exact")
plt.legend(loc="best")
plt.show()








