'''
Created on 1 dec. 2016

@author: EB
'''

import numpy as np
import matplotlib.pyplot as plt
from SQB import *


N = 1000


#===============================================================================
# Cas 1 : volatilite tres faible, temps long, vitesse normale
# On attend une convergence vers la valeur cible
#===============================================================================
T = 10.0
a = 1.5
x = 5.0
t = np.linspace(0,T,N+1)
plt.figure()
plt.plot(t,BESQ_array(a, x, T, N),'b')
plt.grid()
plt.show()
