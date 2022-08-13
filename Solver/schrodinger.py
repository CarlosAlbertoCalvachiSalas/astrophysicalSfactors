import numpy as np
from scipy.integrate import odeint, solve_ivp, ode

e = np.sqrt(1.4399764) # sqrt(MeV . fm)
h = 4.135667696E-21 # MeV. s
hbar = h/(2*np.pi)
c = 299792458

def amuToMeVc_2(m):
	return m*931.49410242 

def amuToMeVfm_2s2(m):
	return (amuToMeVc_2(m)/(c**2))*(1E-30)

def schodingerEquationMatrix(U, E, m, x):
	return np.array([[0, 1],[(2*m*(U(x)-E))/(hbar**2), 0]])

class SchrodingerEquation():
	
	def __init__(self, m, U, radial = False, tIndependent = True):

		self.m = m
		self.U = U
		self.radial = radial
		self.tIndependent = tIndependent

	def solve(self, E, xmin = 0, xmax = 30, Nsteps = 10000):

		M = lambda x: schodingerEquationMatrix(self.U, E, self.m, x)

		diffOperator = lambda x, psi: M(x)@np.array(psi)

		self.k = np.sqrt((2*self.m*E)/(hbar**2))
		self.xLin = np.linspace(xmin, xmax, Nsteps)

		self.T = (2*np.pi)/self.k

		
		return solve_ivp(diffOperator, (0, np.inf),  [1, -1*self.k], method = 'LSODA')
