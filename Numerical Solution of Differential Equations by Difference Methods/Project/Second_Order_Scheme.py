import numpy as np
import Initial as init
import Scheme as sch

# The scheme for second order absorbing boundary conditions

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
            # Special interior points scheme
            U_1[m][n] = dt*init.g(m*dx,n*dy) + a*U_0[m][n] + 0.5*a_x*(U_0[m+1][n] + U_0[m-1][n]) + 0.5*a_y*(U_0[m][n+1]+U_0[m][n-1])

    # Special second order absorbing boundary conditions on all sides:
    for n in range(1,N+1):
        U_1[0][n]=0.5*(2-a_y)*U_0[0][n]-dt*(0.75*np.sqrt(a_x)-1)*init.g(0,n*dy)-0.25*dt*np.sqrt(a_y)*(init.g(2*dx,n*dy)-4*init.g(dx,n*dy))+0.25*a_y*(U_0[0][n+1]+U_0[0][n-1])
        U_1[M+1][n]=0.5*(2-a_y)*U_0[M+1][n]-dt*(0.75*np.sqrt(a_x)-1)*init.g(alpha,n*dy)-0.25*dt*np.sqrt(a_y)*(init.g((M-1)*dx,n*dy)-4*init.g(M*dx,n*dy))+0.25*a_y*(U_0[M+1][n+1]+U_0[M+1][n-1])
    for m in range(1,M+1):
        U_1[m][0] = 0.5*(2-a_x)*U_0[m][0]-dt*(0.75*np.sqrt(a_y)-1)*init.g(m*dx,0)-0.25*dt*np.sqrt(a_x)*(init.g(m*dx,2*dy)-4*init.g(m*dx,dy))+0.25*a_x*(U_0[m+1][0]+U_0[m-1][0])
        U_1[m][N+1] = 0.5*(2-a_x)*U_0[m][N+1]-dt*(0.75*np.sqrt(a_y)-1)*init.g(m*dx,beta)-0.25*dt*np.sqrt(a_x)*(init.g(m*dx,(N-1)*dy)-4*init.g(m*dx,N*dy))+0.25*a_x*(U_0[m+1][N+1]+U_0[m-1][N+1])
    U_1[0][0]=U_0[0][0]-2*dt*(0.25*np.sqrt(a_x)+0.25*np.sqrt(a_y)-0.5)*init.g(0,0)+(dt*np.sqrt(a_y)/6)*(-1*init.g(2*dx,0)+4*init.g(dx,0)) + (dt*np.sqrt(a_x)/6)*(-1*init.g(0,2*dy)+4*init.g(0,dy))
    U_1[0][N+1]=U_0[0][N+1]-2*dt*(0.25*np.sqrt(a_x)+0.25*np.sqrt(a_y)-0.5)*init.g(0,beta)+(dt*np.sqrt(a_y)/6)*(-1*init.g(2*dx,beta)+4*init.g(dx,beta))+(dt*np.sqrt(a_x)/6)*(-1*init.g(0,(N-1)*dy)+4*init.g(0,N*dy))
    U_1[M+1][0]=U_0[M+1][0]-2*dt*(0.25*np.sqrt(a_x)+0.25*np.sqrt(a_y)-0.5)*init.g(alpha,0)+(dt*np.sqrt(a_y)/6)*(-1*init.g((M-1)*dx,0)+4*init.g(M*dx,0))+(dt*np.sqrt(a_x)/6)*(-1*init.g(alpha,2*dy)+4*init.g(alpha,dy))
    U_1[M+1][N+1]=U_0[M+1][N+1]-2*dt*(0.25*np.sqrt(a_x)+0.25*np.sqrt(a_y)-0.5)*init.g(alpha,beta)+(dt*np.sqrt(a_y)/6)*(-1*init.g((M-1)*dx,beta)+4*init.g(M*dx,beta))+(dt*np.sqrt(a_x)/6)*(-1*init.g(alpha,(N-1)*dy)+4*init.g(alpha,N*dy))
    return U_1

def scheme(U_0,U_1,M,N,K,alpha,beta,T):
    dt = T / K
    dx = alpha / (M + 1)
    dy = beta / (N + 1)

    a_x = (dt / dx) ** 2
    a_y = (dt / dy) ** 2

    A_x, A_y = sch.A_matrices(M, N, K, alpha, beta, T)

    # Constants:
    x_const = 1/(0.75*np.sqrt(a_x)+1)
    y_const = 1/(0.75*np.sqrt(a_y)+1)
    const = 1/(1.5+0.75*np.sqrt(a_x)+0.75*np.sqrt(a_y))

    t = dt
    while t < T:
        t += dt
        D = sch.D_matrix(U_1,M,N,K,alpha,beta,T)
        U = np.zeros((M+2,N+2))

        # Matrix equation for the interior points:
        U[1:M+1][:,1:N+1] = -1*U_0[1:M+1][:,1:N+1] + np.matmul(A_x,U_1[1:M+1][:,1:N+1]) + np.matmul(U_1[1:M+1][:,1:N+1],A_y) + D

        # Second order absorbing boundary conditions on all walls:
        for n in range(1,N+1):
            U[0][n] = x_const*((2-a_y)*U_1[0][n]+(0.75*np.sqrt(a_x)-1)*U_0[0][n]+0.25*np.sqrt(a_x)*(-1*U[2][n]+4*U[1][n]+U_0[2][n]-4*U_0[1][n])+0.5*a_y*(U_1[0][n+1]+U_1[0][n-1]))
            U[M+1][n] = x_const*((2-a_y)*U_1[M+1][n]+(0.75*np.sqrt(a_x)-1)*U_0[M+1][n]+0.25*np.sqrt(a_x)*(-1*U[M-1][n]+4*U[M][n]+U_0[M-1][n]-4*U_0[M][n])+0.5*a_y*(U_1[M+1][n+1]+U_1[M+1][n-1]))
        #for m in range(1,M+1):
           # U[m][0] = y_const*((2-a_x)*U_1[m][0]+(0.75*np.sqrt(a_y)-1)*U_0[m][0]+0.25*np.sqrt(a_y)*(-1*U[m][2]+4*U[m][1]+U_0[m][2]-4*U_0[m][1])+0.5*a_y*(U_1[m+1][0]+U_1[m-1][0]))
            #U[m][N+1] = y_const*((2-a_x)*U_1[m][N+1]+(0.75*np.sqrt(a_y)-1)*U_0[m][N+1]+0.25*np.sqrt(a_y)*(-1*U[m][N-1]+4*U[m][N]+U_0[m][N-1]-4*U_0[m][N])+0.5*a_y*(U_1[m+1][N+1]+U_1[m-1][N+1]))
        U[0][0] = const*(3*U_1[0][0]+(0.75*np.sqrt(a_x)+0.75*np.sqrt(a_y)-1.5)*U_0[0][0]+0.25*np.sqrt(a_y)*(-1*U[2][0]+4*U[1][0]+U_0[2][0]-4*U_0[1][0])+0.25*np.sqrt(a_x)*(-1*U[0][2]+4*U[0][1]+U_0[0][2]-4*U_0[0][1]))
        U[0][N+1] = const*(3*U_1[0][N+1]+(0.75*np.sqrt(a_x)+0.75*np.sqrt(a_y)-1.5)*U_0[0][N+1]+0.25*np.sqrt(a_y)*(-1*U[2][N+1]+4*U[1][N+1]+U_0[2][N+1]-4*U_0[1][N+1])+0.25*np.sqrt(a_x)*(-1*U[0][N-1]+4*U[0][N]+U_0[0][N-1]-4*U_0[0][N]))
        U[M+1][0] = const*(3*U_1[M+1][0]+(0.75*np.sqrt(a_x)+0.75*np.sqrt(a_y)-1.5)*U_0[M+1][0]+0.25*np.sqrt(a_y)*(-1*U[M-1][0]+4*U[M][0]+U_0[M-1][0]-4*U_0[M][0])+0.25*np.sqrt(a_x)*(-1*U[M+1][2]+4*U[M+1][1]+U_0[M+1][2]-4*U_0[M+1][1]))
        U[M+1][N+1] = const*(3*U_1[M+1][N+1]+(0.75*np.sqrt(a_x)+0.75*np.sqrt(a_y)-1.5)*U_0[M+1][N+1]+0.25*np.sqrt(a_y)*(-1*U[M-1][N+1]+4*U[M][N+1]+U_0[M-1][N+1]-4*U_0[M][N+1])+0.25*np.sqrt(a_x)*(-1*U[M+1][N-1]+4*U[M+1][N]+U_0[M+1][N-1]-4*U_0[M+1][N]))

        U_0 = U_1
        U_1 = U

    return U