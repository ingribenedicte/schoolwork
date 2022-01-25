import Initial as init
import First_Order_Scheme as scheme_one
import Second_Order_Scheme as scheme_two
import Plots as plots

def main():
    M = 100
    N = 100
    K = 1000

    alpha = 10
    beta = 10
    T = 20
    U_0 = init.initial(M,N,alpha,beta)
    U_1 = scheme_one.first_step(U_0,M,N,K,alpha,beta,T)
    U = scheme_one.scheme(U_0,U_1,M,N,K,alpha,beta,T)

    plots.plot_wave(U,M,N,alpha,beta,T)


main()

