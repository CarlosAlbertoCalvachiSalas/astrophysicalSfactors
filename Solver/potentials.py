import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sfactors.Utils.constants import *
from scipy.integrate import quad

def coulomb(r, Z1, Z2, R = 1):
	if(r < R):
		return ((Z1*Z2*(e**2))/(2*R)) * (3 - ((r/R)**2) )
	else:
		return (Z1*Z2*(e**2))/r

def woodsSaxon(r, V0, R, a):
	return -V0*(1/(1 + np.exp((r - R)/a)))

def expSingle(r, V0, alpha):
	return -V0*np.exp(-alpha*r)

def expDouble(r, V0, alpha0, V1, alpha1):
	return -(V0*np.exp(alpha0*r) + V1*np.exp(alpha1*r))**(-1)

def gaussian(r, V0, alpha):
	return -V0*np.exp(-alpha * (r**2))

def cluster(r, V0, alpha0, V1, alpha1):
	return gaussian(r, V0, alpha0) + V1 * np.exp(-alpha*r)

#def yakovlev(r, Z1, Z2, Rc, Rc1, Ec, beta):

#	U = (Z1*Z2*(e**2)/r)*(r >= Rc1) + (Ec * (1 - beta*(((r - Rc)**2)/(Rc**2))))*(r < Rc1)

#	return U * (U > 0)

#cpotential = coulombPotential(np.linspace(1e-5, 10, 100000), Z1 = 6, Z2 = 6)

r = np.linspace(0, 20, 10000)

#plt.plot(r, np.vectorize(coulomb)(r, Z1 = 6, Z2 = 6), color = 'green')
#plt.plot(r, np.vectorize(coulomb)(r, Z1 = 6, Z2 = 6) + woodsSaxon(r, V0 = 50, R = 10, a = 0.1), color = 'blue')
#plt.plot(r,np.vectorize(coulomb)(r, Z1 = 6, Z2 = 6) + woodsSaxon(r, V0 = 100, R = 5, a = 1), color = 'blue')
#plt.plot(r,np.vectorize(coulomb)(r, Z1 = 6, Z2 = 6), color = 'red')
#plt.plot(r,woodsSaxon(r, V0 = 100, R = 5, a = 1), color = 'green')

#plt.show()
#coulomb(1, 1, 1)
#coulomb(1, *(), **{'Z1': 3, 'Z2': 4})
#coulombPotential(1, *(1,2), **{'Z1': 3, 'Z2': 4})

#print(quad(coulomb, 4, 10, args = (2, 2)))