import numpy as np
import random
import matplotlib.pyplot as plt

# function to return a matrix and vector based on x
def make(x,n):
    cnt = 0
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            A[i][j] = x[cnt]
            A[j][i] = x[cnt]
            cnt += 1
    vec = np.zeros(n)
    for k in range(n):
        vec[k] = x[cnt]
        cnt += 1
    return A, vec

# objective function and gradient for model 1:
def r1(A,c,z,w,m,n):
    r = 0
    if w == 1:
        r = max(np.dot(np.transpose(z-c),np.matmul(A,(z-c)))-1,0)
    if w == -1:
        r = max(1 - np.dot(np.transpose(z-c), np.matmul(A, (z - c))), 0)
    return r

def dr1(A,c,z,w,m,n):
    k = int((n+1)*n/2+n)
    dr = np.zeros(k)
    cnt = 0
    for i in range(n):
        dr[cnt] = (z[i]-c[i])**2
        cnt += 1
        for j in range(i+1,n):
            dr[cnt] = 2*(z[i]-c[i])*(z[j]-c[j])
            cnt += 1
    for l in range(n):
        for k in range(n):
            dr[cnt] += -2*(z[k]-c[k])*(A[k][l])
        cnt += 1
    dr = w*dr
    return dr

def f_1(x,z,w,m,n):
    A,c = make(x,n)
    f = 0
    for i in range(m):
        f += r1(A,c,z[i],w[i],m,n)**2
    return f

def grad_f1(x,z,w,m,n):
    A,c = make(x,n)
    grad_f = 0
    for i in range(m):
        grad_f += 2*r1(A,c,z[i],w[i],m,n)*dr1(A,c,z[i],w[i],m,n)
    return grad_f

# objective function and gradient for model 2:
def r2(A,b,z,w,m,n):
    r = 0
    if w == 1:
        r = max(np.dot(z,np.matmul(A,z)) + np.dot(b,z)-1,0)
    if w == -1:
        r = max(1 - np.dot(z,np.matmul(A,z)) - np.dot(b,z), 0)
    return r

def dr2(A,b,z,w,m,n):
    length = int((n*(n+1))/2 + n)
    dr = np.zeros(length)
    cnt = 0
    for i in range(n):
        for j in range(i,n):
            if i != j:
                dr[cnt] = 2*z[i]*z[j]
            else:
                dr[cnt] = z[i]**2
            cnt += 1
    for l in range(n):
        dr[cnt] = z[l]
        cnt += 1
    dr = w*dr
    return dr

def f_2(x,z,w,m,n):
    A,b = make(x,n)
    f = 0
    for i in range(m):
        f += r2(A,b,z[i],w[i],m,n)**2
    return f

def grad_f2(x,z,w,m,n):
    A,b = make(x,n)
    grad_f = 0
    for i in range(m):
        grad_f += 2*r2(A,b,z[i],w[i],m,n)*dr2(A,b,z[i],w[i],m,n)
    return grad_f

#Question 4, function to verify the implementation:
def test_implementation():
    n = 3
    m = 10
    k = int(n*(n+1)/2 + n)

    # create some random data:
    x = np.random.randn(k)
    p = np.random.randn(k)
    z = np.random.randn(m,n)
    w = np.zeros(m)
    for i in range(n):
        y = random.randint(0,9)
        if y%2 == 0:
            w[i] = 1
        else:
            w[i] = -1

    f1 = f_1(x, z, w, m, n)
    df1 = grad_f1(x, z, w, m, n)
    f2 = f_2(x, z, w, m, n)
    df2 = grad_f2(x, z, w, m, n)

    g1 = df1.dot(p)
    g2 = df2.dot(p)

    print("Model 1:")
    for ep in 10.0**np.arange(-1,-13,-1):
        g_app1 = (f_1(x + ep*p,z,w,m,n)-f1)/ep
        error1 = abs(g_app1-g1)/abs(g1)
        print('ep = %e, error = %e' % (ep,error1))

    print("Model 2:")
    for ep in 10.0**np.arange(-1,-13,-1):
        g_app2 = (f_2(x + ep*p,z,w,m,n)-f2)/ep
        error2 = abs(g_app2-g2)/abs(g2)
        print('ep = %e, error = %e' % (ep,error2))

# line search method which satisfies strong Wolfe:
def LS(model,x,p,z,w,m,n,f,df):
    c_2 = 0.6
    c_1 = 0.5
    a = 1
    a_min = 0
    a_max = np.inf
    cont = True
    while cont:
        if model == "1":
            if f_1(x+a*p,z,w,m,n) > f + c_1*a*np.dot(df,p):
                a_max = a
                a = (a_max + a_min)/2
            elif abs(np.dot(grad_f1(x+a*p,z,w,m,n),p)) > c_2*abs(np.dot(df,p)):
                a_min = a
                if np.isinf(a_max):
                    a = 2*a
                else:
                    a = (a_max + a_min)/2
            else:
                cont = False
        else:
            if f_2(x+a*p,z,w,m,n) > f + c_1*a*np.dot(df,p):
                a_max = a
                a = (a_max + a_min)/2
            elif abs(np.dot(grad_f2(x+a*p,z,w,m,n),p)) > c_2*abs(np.dot(df,p)):
                a_min = a
                if np.isinf(a_max):
                    a = 2*a
                else:
                    a = (a_max + a_min)/2
            else:
                cont = False
    return a

# Fletcher Reeves method:
def FR(z, w, n, m, model, x, tol):
    if model == "1":
        df = grad_f1(x,z,w,m,n)
        f = f_1(x,z,w,m,n)
    else:
        df = grad_f2(x, z, w, m, n)
        f = f_2(x,z,w,m,n)
    p = -1*df
    cnt = 0
    while np.linalg.norm(df,ord=np.inf) > tol and f != 0:
        a = LS(model,x,p,z,w,m,n,f,df)
        x_next = x + a*p
        if model == "1":
            df_next = grad_f1(x_next,z,w,m,n)
            f = f_1(x_next,z,w,m,n)
        else:
            df_next = grad_f2(x_next,z,w,m,n)
            f = f_2(x_next,z,w,m,n)
        B = np.dot(np.transpose(df_next),df_next) / np.dot(np.transpose(df),df)
        p_next = -1*df_next + B*p

        x = x_next
        p = p_next
        df = df_next
        cnt = cnt + 1
    return np.linalg.norm(df,ord=np.inf),x, cnt

# Steepest descent method:
def SD(z, w, n, m, model, x, tol):
    if model == "1":
        df = grad_f1(x,z,w,m,n)
        f = f_1(x,z,w,m,n)
    else:
        df = grad_f2(x,z,w,m,n)
        f = f_2(x,z,w,m,n)
    cnt = 0
    while np.linalg.norm(df,ord=np.inf) > tol:
        a = LS(model,x,-1*df,z,w,m,n,f,df)
        x = x - a*df
        if model == "1":
            df = grad_f1(x, z, w, m, n)
            f = f_1(x,z,w,m,n)
        else:
            df = grad_f2(x, z, w, m, n)
            f = f_2(x,z,w,m,n)
        cnt = cnt + 1
    return np.linalg.norm(df,ord=np.inf), x, cnt

# Testing the optimization algorithms:
def Q5():
    m = 100
    n = 20
    k = int((n+1)*n/2+n)
    tol = 10**(-5)

    x_0 = np.ones(k)
    z = np.random.randn(m,n)
    w = np.zeros(m)
    for i in range(n):
        y = random.randint(0,9)
        if y%2 == 0:
            w[i] = 1
        else:
            w[i] = -1

    FR_1,x_1,iterations_FR1 = FR(z,w,n,m,"1",x_0,tol)
    print("FR-method, model 1")
    print("error: ",FR_1,"        number of iterations: ",iterations_FR1)
    SD_1, x_11, iterations_SD1 = SD(z,w,n,m,"1",x_0,tol)
    print("SD-method, model 1")
    print("error: ",SD_1,"        number of iterations: ",iterations_SD1)

    FR_2,x_2,iterations_2 = FR(z,w,n,m,"2",x_0,tol)
    print("FR-method, model 2")
    print("error: ", FR_2,"       number of iterations: ",iterations_2)
    SD_2, x_11, iterations_SD2 = SD(z,w,n,m,"2",x_0,tol)
    print("SD-method, model 2")
    print("error: ",SD_2,"        number of iterations: ",iterations_SD1)

Q5()

# Approximation to 2-dimensional rectangle set:
def make_ellipse(A,vec,model,area):
    x = np.arange(-area, area, 0.1)
    y = np.arange(-area, area, 0.1)
    X, Y = np.meshgrid(x, y)
    if model == "1":
        Z = A[0][0]*(X-vec[0])**2 + 2*A[0][1]*(X-vec[0])*(Y-vec[1]) + A[1][1]*(Y-vec[1])**2
    else:
        Z = A[0][0]*X**2 + 2*A[0][1]*X*Y + A[1][1]*Y**2 + vec[0]*X + vec[1]*Y
    return X,Y,Z

def plot_z(z,w,m):
    for i in range(m):
        if w[i]<0:
            col='r'
        else:
            col='b'
        plt.plot(z[i][0], z[i][1], 'x', color=col)
        
def approximate_rectangle_2D():
    m = 200
    n = 2
    k = int((n+1)*n/2+n)
    tol = 10**(-5)
    area = 1.0
    rectangle = [-area/np.random.uniform(1,2),area/np.random.uniform(1,2),-area/np.random.uniform(1,2),area/np.random.uniform(1,2)]
    x_0 = np.ones(k)
    z = np.zeros((m,n))
    w = np.zeros(m)
    for i in range(m):
        z[i] = np.random.uniform(-area,area,n)
        x = z[i][0]
        y = z[i][1]
        if rectangle[0] < x < rectangle[1] and rectangle[2] < y < rectangle[3]:
            w[i] = 1
        else: 
            w[i] = -1

    df_FR1,x_FR1,cnt_FR1 = FR(z,w,n,m,"1",x_0,tol)
    A,c = make(x_FR1,n)
    X,Y,Z = make_ellipse(A,c,"1",area)
    plt.figure()
    plt.title("The FR-method, model 1")
    set = plt.contour(X,Y,Z,[1])
    plt.clabel(set)
    plot_z(z,w,m)
    plt.show()
    df_FR2,x_FR2,cnt_FR2 = FR(z,w,n,m,"2",x_0,tol)
    A,b = make(x_FR2,n)
    X,Y,Z = make_ellipse(A,b,"2",area)
    plt.figure()
    plt.title("The FR-method, model 2")
    set = plt.contour(X,Y,Z,[1])
    plt.clabel(set)
    plot_z(z,w,m)
    plt.show()
    df_SD1,x_SD1,cnt_SD1 = SD(z,w,n,m,"1",x_0,tol)
    A,c = make(x_SD1,n)
    X,Y,Z = make_ellipse(A,c,"1",area)
    plt.figure()
    plt.title("The SD-method, model 1")
    set = plt.contour(X,Y,Z,[1])
    plt.clabel(set)
    plot_z(z,w,m)
    plt.show()
    df_SD2,x_SD2,cnt_SD2 = SD(z,w,n,m,"2",x_0,tol)
    A,b = make(x_SD2,n)
    X,Y,Z = make_ellipse(A,b,"2",area)
    plt.figure()
    plt.title("The SD-method, model 2")
    set = plt.contour(X,Y,Z,[1])
    plt.clabel(set)
    plot_z(z,w,m)
    plt.show()