import pylab as py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from WishartDim2_O2 import *

#parametres
T = 1. #horizon
a = 1.1
x = np.array([[1,0],[0,1]])
N = 1 #nombre initial de pas de discretisation
M = int(1.e2) #nombre de simulations independantes
P = 6 #nombre d'iterations sur la valeur du pas
reals = np.zeros(M) #pour mettre les tirages des Wishart

#matrice z fixee, parametres u et v    
z = np.array([[0.04,0.02],[0.02,0.04]])
u = 0.2
v = 0.5 

#Fonction caracteristique cas a = u*I2 et b = v*I2
def LaplaceTransform(z,t,a,x,u,v):
    q = np.array([[u**2/(2*v)*(np.exp(2*v*t)-1),0],
                  [0,u**2/(2*v)*(np.exp(2*v*t)-1)]])
    m = np.array([[np.exp(t*v),0],
                  [0,np.exp(t*v)]])
    e = np.eye(2)-2*py.dot(q,z)
    w = py.inv(e)
    w = np.dot(w,np.dot(m,np.dot(x,m)))
    y = py.dot(z,w)
    den = py.det(e)**(a/2.0)
    num = np.exp(py.trace(y))
    return num/den
    
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

err_eul = np.zeros(P) #vecteur des erreurs faibles du schema d'Euler
largeur_IC_eul = np.zeros(P) #demi-largeur des intervalles de confiance a 95%
#contr_err_eul = np.zeros(P) #vecteur des erreurs faibles Euler avec variable de controle
contr_largeur_IC_eul = np.zeros(P) #demi-largeur des intervalles de confiance 95% avec var. de controle
Npas = np.zeros(P) # vecteur des nombres de pas

lap_exact = LaplaceTransform_e1(z,T,a,x)

for i in range(P): #boucle qui incremente le nombre de pas
  
    sum_laplace = 0
    sum_laplace_sqr = 0

    for k in range(M): #boucle sur les pas de temps
  
        wis = WIS2_O2_array_e1(a,x,T,N)
        lap_emp = expoTrace(z,wis[N])
        sum_laplace += lap_emp
        sum_laplace_sqr += lap_emp**2

  
    err_eul[i] = sum_laplace/M - lap_exact
    largeur_IC_eul[i] = 1.96 * np.sqrt((sum_laplace_sqr/M - (sum_laplace/M)**2)/M)

    Npas[i] = N
    N = N*2 #multiplication du nombre N de pas par 2

print("Erreurs faibles Euler: "); print(err_eul)
print("demi-largeur IC 95%: "); print(largeur_IC_eul);  print("\n")

#representation graphique de (h=T/N, contr_err_eul)
plt.clf()
plt.plot(T/Npas, err_eul, color="r", label="erreur faible")
plt.plot(T/Npas, err_eul - largeur_IC_eul, color="b", label="IC 95%")
plt.plot(T/Npas, err_eul + largeur_IC_eul, color="b")
plt.xlabel('Pas de discretisation')
plt.legend(loc="best")
plt.show()