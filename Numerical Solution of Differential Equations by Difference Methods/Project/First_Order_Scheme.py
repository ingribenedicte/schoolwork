import numpy as np
import Initial as init
import Scheme as sch

# The scheme for first order absorbing boundary conditions

# Computes the first time step:
def first_step(U_0,M,N,K,alpha,beta,T):
    dt = T/K
    dx = alpha/(M+1)
    dy = beta/(N+1)

    a_x = (dt/dx)**2
    a_y = (dt/dy)**2
    a = 1 - a_x - a_y

    U_1 = np.zeros((M+2,N+2))
    for m in range(1,M+1):
        for n in range(1,N+1):
            # Special interior points scheme:
            U_1[m][n] = dt*init.g(m*dx,n*dy) + a*U_0[m][n] + 0.5*a_x*(U_0[m+1][n] + U_0[m-1][n]) + 0.5*a_y*(U_0[m][n+1]+U_0[m][n-1])

    # Special first order absorbing boundary conditions on all sides:
    for n in range(1,N+1):
        U_1[0][n] = a*U_0[0][n]-(np.sqrt(a_x)-1)*dt*init.g(0,n*dy)+a_x*U_0[1][n]+0.5*a_y*(U_0[0][n+1]+U_0[0][n-1])
        U_1[M+1][n] = a*U_0[M+1][n]-(np.sqrt(a_x)-1)*dt*init.g(alpha,n*dy)+a_x*U_0[M][n]+0.5*a_y*(U_0[M+1][n+1]+U_0[M+1][n-1])
    for m in range(1,M+1):
        U_1[m][0] = a*U_0[m][0]-(np.sqrt(a_y)-1)*dt*init.g(m*dx,0)+a_y*U_0[m][1]+0.5*a_x*(U_0[m+1][0]+U_0[m-1][0])
        U_1[m][N+1] = a*U_0[m][N+1]-(np.sqrt(a_y)-1)*dt*init.g(m*dx,beta)+a_y*U_0[m][N]+0.5*a_x*(U_0[m+1][N+1]+U_0[m-1][N+1])
    U_1[0][0] = (1-np.sqrt(a_x)-np.sqrt(a_y))*dt*init.g(0,0)+a*U_0[0][0]+a_x*U_0[1][0]+a_y*U_0[0][1]
    U_1[0][N+1] = (1-np.sqrt(a_x)-np.sqrt(a_y))*dt*init.g(0,beta)+a*U_0[0][N+1]+a_x*U_0[1][N+1]+a_y*U_0[0][N]
    U_1[M+1][0] = (1-np.sqrt(a_x)-np.sqrt(a_y))*dt*init.g(alpha,0)+a*U_0[M+1][0]+a_x*U_0[M][0]+a_y*U_0[M+1][1]
    U_1[M+1][N+1] = (1-np.sqrt(a_x)-np.sqrt(a_y))*dt*init.g(alpha,beta)+a*U_0[M+1][N+1]+a_x*U_0[M][N+1]+a_y*U_0[M+1][N]
    return U_1


def scheme(U_0,U_1,M,N,K,alpha,beta,T):
    dt = T / K
    dx = alpha / (M + 1)
    dy = beta / (N + 1)

    a_x = (dt / dx) ** 2
    a_y = (dt / dy) ** 2
    a = 1 - a_x - a_y

    A_x,A_y = sch.A_matrices(M,N,K,alpha,beta,T)

    t = dt
    while t < T:
        t += dt
        D = sch.D_matrix(U_1,M,N,K,alpha,beta,T)
        U = np.zeros((M+2,N+2))

        # Matrix equation for the interior points:
        U[1:M+1][:,1:N+1] = -1*U_0[1:M+1][:,1:N+1] + np.matmul(A_x,U_1[1:M+1][:,1:N+1]) + np.matmul(U_1[1:M+1][:,1:N+1],A_y) + D

        # First order absorbing boundary conditions on all walls:
        for n in range(1,N+1):
            U[0][n] = (1/(np.sqrt(a_x)+1))*(2*a*U_1[0][n]+2*a_x*U_1[1][n]+(np.sqrt(a_x)-1)*U_0[0][n]+a_y*(U_1[0][n+1]+U_1[0][n-1]))
            U[M+1][n] = (1/(np.sqrt(a_x)+1))*(2*a*U_1[M+1][n]+2*a_x*U_1[M][n]+(np.sqrt(a_x)-1)*U_0[M+1][n]+a_y*(U_1[M+1][n+1]+U_1[M+1][n-1]))
        for m in range(1,M+1):
            U[m][0] = (1/(np.sqrt(a_y)+1))*(2*a*U_1[m][0]+2*a_y*U_1[m][1]+(np.sqrt(a_y)-1)*U_0[m][0]+a_y*(U_1[m+1][0]+U_1[m-1][0]))
            U[m][N+1] = (1/(np.sqrt(a_y)+1))*(2*a*U_1[m][N+1]+2*a_y*U_1[m][N]+(np.sqrt(a_y)-1)*U_0[m][N+1]+a_y*(U_1[m+1][N+1]+U_1[m-1][N+1]))
        U[0][0] = (1/(1+np.sqrt(a_x)+np.sqrt(a_y)))*((-1+np.sqrt(a_x)+np.sqrt(a_y))*U_0[0][0]+2*a*U_1[0][0]+2*a_x*U_1[1][0]+2*a_y*U_1[0][1])
        U[0][N+1] = (1/(1+np.sqrt(a_x)+np.sqrt(a_y)))*((-1+np.sqrt(a_x)+np.sqrt(a_y))*U_0[0][N+1]+2*a*U_1[0][N+1]+2*a_x*U_1[1][N+1]+2*a_y*U_1[0][N])
        U[M+1][0] = (1/(1+np.sqrt(a_x)+np.sqrt(a_y)))*((-1+np.sqrt(a_x)+np.sqrt(a_y))*U_0[M+1][0]+2*a*U_1[M+1][0]+2*a_x*U_1[M][0]+2*a_y*U_1[M+1][1])
        U[M+1][N+1] = (1/(1+np.sqrt(a_x)+np.sqrt(a_y)))*((-1+np.sqrt(a_x)+np.sqrt(a_y))*U_0[M+1][N+1]+2*a*U_1[M+1][N+1]+2*a_x*U_1[M][N+1]+2*a_y*U_1[M+1][N])

        U_0 = U_1
        U_1 = U
    return U