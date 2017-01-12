import numpy as np

# def SQB(a,n,T,x):
#     t = np.linspace(0,T,n+1)
#     mu = a/2.0 - 1
#     r = np.zeros(n+1)
#     r[0] = x
#     for i in np.arange(n):
#         Y = np.random.poisson(r[i]/(2*(t[i+1]-t[i])))
#         r[i+1] = np.random.gamma(Y+mu+1,1.0/(2*(t[i+1]-t[i])))
#     return r

#Simulation d'un processus de Bessel carre BESQ
#verifiant l'EDS suivante : dXt = n*dt + 2*sqrt(Xt)dWt
#ou n est un reel quelconque
#instant initial : x

def BESQ_array(a,x,T,N):
    t = np.linspace(0,T,N+1)
    r = np.zeros(N+1)
    r[0] = x
    if (a>1):
        for i in np.arange(N):
            lam = r[i]/(t[i+1]-t[i])
            c = (t[i+1]-t[i])
            Z = np.random.randn()
            X = np.random.chisquare(a-1)
            r[i+1] = c*((Z+np.sqrt(lam))**2+X)
    else:
        for i in np.arange(N):
            lam = r[i]/(t[i+1]-t[i])
            c = (t[i+1]-t[i])
            P = np.random.poisson(lam/2)
            X = np.random.chisquare(a+2*P)
            r[i+1] = c*X
    return r
    
#BESQ pour passer de t1 a t2
def BESQ(a,x,t_1,t_2):
    if (a>1):
        lam = x/(t_2-t_1)
        c = (t_2-t_1)
        Z = np.random.randn()
        X = np.random.chisquare(a-1)
        return c*((Z+np.sqrt(lam))**2+X)
    else:
        lam = x/(t_2-t_1)
        c = (t_2-t_1)
        P = np.random.poisson(lam/2)
        X = np.random.chisquare(a+2*P)
        return c*X
    

