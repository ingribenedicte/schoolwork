import numpy as np

# function 1
f1 = lambda x: x**2 - 4*x + 3
df1 = lambda x: 2*x - 4 # function 1 derived
ddf1 = lambda x: 2 # function 1 double derived
root1_1 = 1
root1_2 = 3

# function 2
f2 = lambda x: np.exp(2*x) - 6*np.exp(x) + 8
df2 = lambda x: 2*np.exp(2*x) - 6*np.exp(x) # function 2 derived
ddf2 = lambda x: 4*np.exp(2*x) - 6*np.exp(x) # function 2 double derived
root2_1 = np.log(2)
root2_2 = 2*np.log(2)

# function 3
f3 = lambda x: x**2 * np.exp(-x**2)
df3 = lambda x: 2 * x * np.exp(-x**2) - 2* x**3 * np.exp(-x**2) # function 3 derived
ddf3 = lambda x: 4* x**4 * np.exp(-x**2) - 10 * x**2 - np.exp(-x**2) + 2 * np.exp(-x**2) # function 3 double derived
root3 = 0