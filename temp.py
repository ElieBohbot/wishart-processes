#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 21:00:35 2017

@author: EB
"""
import numpy as np
from WishartDim2 import *
import pylab as py

#Parametres de test
T = 1.0
a = 1.1
x = np.array([[1,0],[0,1]])
M = int(1.e1)

#matrice z fixee, parametres u et v    
z = np.array([[0.04,0.02],[0.02,0.04]])
u = 0.2
v = 0.5 

simul = np.zeros((M,4,2,2))
