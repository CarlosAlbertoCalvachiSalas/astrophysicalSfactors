import numpy as np
from scipy.integrate import quad, nquad
from sfactors.Solver.potentials import gaussian

"""
For details related to parametrization see:

F. Koyuncu and A. Soylu,
Chinese Physics C 42, 054106 (2018). 

"""

def densityM(r, rho0M, omega, betaM):
	return rho0M*(1+omega*(r**2))*np.exp(-betaM*(r**2))

def densityAlpha(r, rho0Alpha, betaAlpha):
	return rho0Alpha*np.exp(-betaAlpha*(r**2))

def sphericalVector(r, theta, phi):
	return np.array([	r*np.sin(theta)*np.cos(phi),
						r*np.sin(theta)*np.sin(phi), 
						r*np.cos(theta)])

def r12(r, theta, phi, r1, theta1, phi1, r2, theta2, phi2):

	R 		= sphericalVector(r, theta, phi)
	r1Vect	= sphericalVector(r1, theta1, phi1)
	r2Vect	= sphericalVector(r2, theta2, phi2)

	r12Vect = R - r1Vect + r2Vect

	return np.linalg.norm(r12Vect)

def doubleFoldingFunction(r, theta, phi, rho0Alpha, rho0M, omega, betaAlpha, betaM, v0, gamma):

	def f(r1, theta1, phi1, r2, theta2, phi2): 

		r12Distance = r12(r, theta, phi, r1, theta1, phi1, r2, theta2, phi2)
		weight = (r1**2) * np.sin(theta1) * (r2**2) * np.sin(theta2)

		return densityM(r, rho0M, omega, betaM)*densityAlpha(r, rho0Alpha, betaAlpha)*gaussian(r12Distance, v0, gamma)*weight

	return np.vectorize(f)

def doubleFolding(r, theta, phi, rho0Alpha, rho0M, omega, betaAlpha, betaM, v0, gamma):

	f = doubleFoldingFunction(r, theta, phi, rho0Alpha, rho0M, omega, betaAlpha, betaM, v0, gamma)

	return nquad(f, [	[0, 1], [0, np.pi], [0, 2*np.pi], 
						[0, 1], [0, np.pi], [0, 2*np.pi] ])

"""
doubleFolding(r = 0, 
			 theta = 0, 
			 phi = 0, 
			 rho0Alpha = 0.4229, 
			 rho0M = 0.1644, 
			 omega = 0.4988, 
			 betaAlpha = 0.7024,
			 betaM = 0.3741,
			  v0 = 122.6225, 
			  gamma = 0.22)
"""

def test(r1, theta1, phi1, r2, theta2, phi2):

	r1vec = np.array([	r1*np.sin(theta1)*np.cos(phi1),
						r1*np.sin(theta1)*np.sin(phi1), 
						r1*np.cos(theta1)])

	r2vec = np.array([	r2*np.sin(theta2)*np.cos(phi2),
						r2*np.sin(theta2)*np.sin(phi2), 
						r2*np.cos(theta2)])

	r12Distance = np.linalg.norm(r1vec - r2vec)

	return np.exp(-r12Distance**2)

nquad(test, [	[0, 1], [0, np.pi], [0, 2*np.pi], 
						[0, 1], [0, np.pi], [0, 2*np.pi] ])
