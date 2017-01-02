'''
Created on 1 dec. 2016

@author: EB
'''

import numpy as np
import matplotlib.pyplot as plt
from CIR import *


n = 8


#===============================================================================
# Cas 1 : volatilite tres faible, temps long, vitesse normale
# On attend une convergence vers la valeur cible
#===============================================================================
T = 20.0
a = 0.2
b = 0.05
sig = 0.001
t = np.linspace(0,T,2**n+1)
plt.figure()
plt.plot(t,CIR(a, b, sig, n, T),'b')
plt.grid()
plt.show()

#===============================================================================
# Cas 2 : volatilite normal, temps intermedaire, vitesse tres forte
# On attend une stabilisation rapide autour de la valeur cible
#===============================================================================
T = 5.0
a = 20
b = 0.05
sig = 0.1
t = np.linspace(0,T,2**n+1)
plt.figure()
plt.plot(t,CIR(a, b, sig, n, T),'b')
plt.grid()
plt.show()
