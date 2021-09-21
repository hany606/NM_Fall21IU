from scipy.integrate import quad
import numpy as np

def integrand(x, i):
    l = [2*x+1, np.cos(x), np.sin(x), np.tan(x), x*x, np.exp(x), np.log10(x+2), np.sqrt(x+1),  x*np.tan(x)*np.tan(x), np.exp(-x*5), np.exp(-x*x), np.sin(x*x*x), np.sin(x*x), x**10]
    return l[i]#x**10#np.pow(x,10)#np.sin(x*x*x)#np.exp(-x*x)

a = -1
b = 1

for i in range(14):
    I = quad(integrand, a, b, args=(i))
    print(I)