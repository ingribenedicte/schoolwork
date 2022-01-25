import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import Error as error

# Surface plot of wave, works for numerical and analytical:
def plot_wave(U, M, N, alpha, beta, T):
    [xx, yy] = np.meshgrid(np.linspace(0, alpha, N + 2), np.linspace(0, beta, M + 2))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(xx, yy, U, rstride=1, cstride=1, alpha=1, linewidth=1, cmap=cm.hot, antialiased=True)
    ax.set_zlim(-1,1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis([0, alpha, 0, beta])
    title = "The wave at t = " + str(T)
    plt.title(title)
    plt.show()

def plot_spatial_error(K, alpha, beta, T):
    steplength_list,error_list = error.spatial_error(K,alpha,beta,T)
    plt.figure()
    plt.title("Spatial error")
    plt.xlabel("Spatial step length")
    plt.ylabel("Infinity norm")
    plt.loglog(steplength_list, error_list)
    plt.show()