import numpy as np
import math
import Scheme as scheme

# Analytical solution:
def u(M,N,alpha,beta,t):
    dx = alpha/(M+1)
    dy = beta/(N+1)
    u = np.zeros((M+2,N+2))
    for m in range(M+2):
        for n in range(N+2):
            u[m][n] = np.sin(np.pi*m*dx)*np.sin(np.pi*n*dy)*np.cos(np.sqrt(2)*np.pi*t)
    return u

# Error in space:
def spatial_error(K,alpha,beta,T):
    steplength_list = np.zeros(5)
    error_list = np.zeros(5)
    for i in range(5):
        M = 2**(i+3)
        N = 2**(i+3)
        dx = alpha/(M+1)
        dy = beta/(N+1)

        u_exact = u(M,N,alpha,beta,T)
        U = scheme.scheme(M,N,K,alpha,beta,T)

        steplength_list[i] = dx
        error = np.subtract(U,u_exact)
        error_list[i] = np.sqrt(dy*dx)*np.linalg.norm(error,2)
    return steplength_list,error_list

# Error in time:
def error_time(alpha,beta,T):
    steplength_list = np.zeros(7)
    error_list = np.zeros(7)
    for i in range(7):
        M = 2**(i+2)
        N = 2**(i+2)
        K = math.ceil(np.sqrt(2)*10*(M+1))

        steplength_list[i] = T/K

        U_0 = init.initial(M,N,alpha,beta)
        U_1 = scheme.first_step(U_0,M,N,K,alpha,beta,T)
        U = scheme_reg(U_0,U_1,M,N,K,alpha,beta,T)
        u_exact = u(M,N,alpha,beta,T)

        error = np.subtract(U,u_exact)
        error_list[i] = np.linalg.norm(error,ord=np.inf)
    return steplength_list,error_list