'''
Created on 1 dec. 2016

@author: EB
'''
import numpy as np

#===============================================================================
# Simulation d'un processus CIR par l'algorithme de Glasserman (loi exacte)
# Grille [0,T] subdivisee en 2**n+1 points 
#===============================================================================
def CIR(a,b,sig,n,T):
    t = np.linspace(0,T,2**n+1)
    d = 4*b*a/(sig**2)
    r = np.zeros(2**n+1)
    r[0] = 0.8*b
    if (d>1):
        for i in np.arange(2**n):
            c = sig**2*(1-np.exp(-a*(t[i+1]-t[i])))/(4*a)
            lam = r[i]*np.exp(-a*(t[i+1]-t[i]))/c
            Z = np.random.randn()
            X = np.random.chisquare(d-1)
            r[i+1] = c*((Z+np.sqrt(lam))**2+X)
    else:
        for i in np.arange(2**n):
            c = sig**2*(1-np.exp(-a*(t[i+1]-t[i])))/(4*a)
            lam = r[i]*np.exp(-a*(t[i+1]-t[i]))/c
            N = np.random.poisson(lam/2)
            X = np.random.chisquare(d+2*N)
            r[i+1] = c*X
    return r
