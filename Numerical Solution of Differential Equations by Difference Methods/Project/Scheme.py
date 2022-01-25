import numpy as np
import Initial as init

# This file holds the explicit leapfrog scheme without absorbing boundary conditions

def first_step(U_0, M, N, K, alpha, beta, T):
    dt = T/K
    dx = alpha/(M+1)
    dy = beta/(N+1)

    a_x = (dt/dx)**2
    a_y = (dt/dy)**2
    a = 1 - a_x - a_y

    U_1 = np.zeros((M+2,N+2))
    for m in range(1,M+1):
        for n in range(1,N+1):
            U_1[m][n] = dt*init.g(m*dx,n*dy) + a*U_0[m][n] + 0.5*a_x*(U_0[m+1][n] + U_0[m-1][n]) + 0.5*a_y*(U_0[m][n+1]+U_0[m][n-1])

    return U_1

def D_matrix(U,M,N,K,alpha,beta,T):
    dt = T/K
    dx = alpha/(M+1)
    dy = beta/(N+1)

    s_x = (dt/dx)**2
    s_y = (dt/dy)**2

    D = np.zeros((M,N))
    D[0][0] = s_x*U[0][1] + s_y*U[1][0]
    D[0][N-1] = s_x*U[0][N] + s_y*U[1][N+1]
    D[M-1][0] = s_x*U[M+1][1] + s_y*U[M][0]
    D[M-1][N-1] = s_x*U[M+1][N] + s_y*U[M][N+1]
    for m in range(1,M-1):
        D[m][0] = s_x*U[m+1][0]
        D[m][N-1] = s_x*U[m+1][N+1]
    for n in range(1,N-1):
        D[0][n] = s_y*U[0][n+1]
        D[M-1][n] = s_y*U[M+1][n+1]
    return D

def A_matrices(M,N,K,alpha,beta,T):
    dt = T/K
    dx = alpha/(M+1)
    dy = beta/(N+1)

    a_x = (dt/dx)**2
    a_y = (dt/dy)**2
    a = 1 - a_x - a_y

    A_x = np.zeros((M,M))
    for m in range(M):
        A_x[m][m] = a
        if m != 0:
            A_x[m][m-1] = a_x
        if m != M-1:
            A_x[m][m+1] = a_x

    A_y = np.zeros((N, N))
    for n in range(N):
        A_y[n][n] = a
        if n != 0:
            A_y[n][n-1] = a_y
        if n != N-1:
            A_y[n][n+1] = a_y

    return A_x, A_y

def scheme(M,N,K,alpha,beta,T):
    dt = T/K

    A_x, A_y = A_matrices(M, N, K, alpha, beta, T)
    U_0 = init.initial(M, N, alpha, beta)
    U_1 = first_step(U_0, M, N, K, alpha, beta, T)

    t = dt
    while t < T:
        t += dt
        D = D_matrix(U_1,M,N,K,alpha,beta,T)
        U = np.zeros((M+2,N+2))

        U[1:M+1][:,1:N+1] = -1*U_0[1:M+1][:,1:N+1] + np.matmul(A_x,U_1[1:M+1][:,1:N+1]) + np.matmul(U_1[1:M+1][:,1:N+1],A_y) + D

        U_0 = U_1
        U_1 = U
    return U