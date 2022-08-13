import numpy as np
from scipy.special import hyp1f1, hyperu, gamma
from scipy.integrate import solve_ivp
from mpmath import coulombf, coulombg

def coulombDiff(eta, l):

	f = lambda rho, psi: np.array([psi[1], -((1 - (2*eta/rho) - ((l*(l+1))/rho**2))*psi[0])])

	return solve_ivp(f, (1e-12, 20),  [6, -1e5], t_eval= np.linspace(1e-12, 20, 200))

def coulombHpos(rho, eta, l):
	C = ((2**l)*np.sqrt(gamma(l + 1 + (1j*eta))*gamma(l + 1 - (1j*eta))))/(gamma(2*l + 2)*np.exp(eta*(np.pi/2)))

	Dpos = ((-2j)**(2*l+1))*(gamma(l + 1 + (1j*eta))/(C*gamma(2*l + 2)))

	#print(hyp1f1(l + 1 + (1j*eta), 1, 1))
	return Dpos * (rho**(l+1)) * np.exp(1j*rho) * hyperu(l + 1 + (1j*eta), 2*l + 2,-(2j*rho))

def coulombHneg(rho, eta, l):
	C = ((2**l)*np.sqrt(gamma(l + 1 + (1j*eta))*gamma(l + 1 - (1j*eta))))/(gamma(2*l + 2)*np.exp(eta*(np.pi/2)))
	
	Dneg = ((2j)**(2*l+1))*(gamma(l + 1 - (1j*eta))/(C*gamma(2*l + 2)))

	return Dneg * (rho**(l+1)) * np.exp(-1j*rho) * hyperu(l + 1 - (1j*eta), 2*l + 2,(2j*rho))

def coulombF(rho, eta, l):
	
	#return np.imag(coulombHpos(rho, eta, l))

	return np.vectorize(coulombf, excluded = [0,1])(l, eta, rho)

def coulombG(rho, eta, l):
	
	return np.vectorize(coulombg, excluded = [0,1])(l, eta, rho)

	#return np.real(coulombHpos(rho, eta, l))



