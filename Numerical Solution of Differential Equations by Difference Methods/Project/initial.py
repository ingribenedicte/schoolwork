import numpy as np

# initial conditions:
def f(x,y):
    return np.exp(-0.5*((x-5)**2+(y-5)**2))
    #return np.sin(np.pi*x)*np.sin(np.pi*y)

def g(x,y):
    return 0

def initial(M,N,alpha,beta):
    dx = alpha/(M+1)
    dy = beta/(N+1)
    U_0 = np.zeros((M+2,N+2))
    for m in range(M+2):
        for n in range(N+2):
            U_0[m][n] = f(m*dx,n*dy)
    return U_0
