'''
Created on 1 dec. 2016

@author: EB
'''

import numpy as np
import matplotlib.pyplot as plt
from SQB import *


n = 8

#===============================================================================
# Cas 1 : volatilite tres faible, temps long, vitesse normale
# On attend une convergence vers la valeur cible
#===============================================================================
T = 5.0
a = 3.0
x = 5.0
t = np.linspace(0,T,2**n+1)
plt.figure()
plt.plot(t,BESQ(a, x, T, n),'b')
plt.grid()
plt.show()
