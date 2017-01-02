#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 18:25:19 2017

@author: EB
"""

import numpy as np
from SQB import *

N = 1000 #toutes les simulations auront N+1 points
a = 1.5
x = np.array([[2,1],[1,1]])

#==============================================================================
# Algorithm 1 : Exact simulation WIS_2(x,a,0,e^1_2;t)
#==============================================================================

def Algo1(x,a):
    X = np.zeros(N+1,(2,2))
    if (x[1][1]==0):
        X[0]=x
        


