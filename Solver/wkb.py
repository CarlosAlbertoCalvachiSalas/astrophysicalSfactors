import numpy as np
from scipy.optimize import fsolve, fmin
from scipy.integrate import quad, quad_vec
from scipy.misc import derivative
import matplotlib.pyplot as plt

from sfactors.nonResonant.yakovlev import *

e = np.sqrt(1.4399764) # sqrt(MeV . fm)
h = 4.135667696E-21 # MeV. s
hbar = h/(2*np.pi)

def amuToMeVc_2(m):
	return m*931.49410242 

def amuToMeVfm_2s2(m):
	return (amuToMeVc_2(m)/(c**2))*(1E-30)


class WKB():

	def __init__(self, U, m, rmin = 0, rmax = 20, Emin = 0, Emax = 30, Ncrit = 10000):	
		self.U = U # U(r)
		self.invU = lambda r: -self.U(r)
		self.vectU = np.vectorize(self.U, otypes = [float])
		self.vectInvU = np.vectorize(self.invU, otypes = [float])
		self.Emin = Emin
		self.Emax = Emax

		self.m = m

		self.rleft = rmin
		self.right = rmax
		self.Ncrit = Ncrit
		self.rLin = np.linspace(self.rleft, self.right, self.Ncrit)

	def graphPotential(self):

		plt.plot(self.rLin, self.vectU(self.rLin), color = 'blue', label = 'plot')
		plt.ylabel('U(MeV)')
		plt.xlabel('r(fm)')
		plt.title('Potential plot')

	def f(self, E):

		K = - (2/hbar)

		self.phi = lambda r: K * np.sqrt((2*self.m)*(self.vectU(r) - E)) 
		return self.phi 

	def phiSinge(self, E, r1, r2):

		return np.array(quad(self.f(E), r1, r2))

	def phi(self, E, r1, r2):

		excludedParams = np.array([])

		if(not hasattr(E, '__len__')):
			excludedParams = np.append([0], excludedParams)

		if(not hasattr(r1, '__len__')):
			excludedParams = np.append([1], excludedParams)

		if(not hasattr(r2, '__len__')):
			excludedParams = np.append([2], excludedParams)

		signatureString = '(),'*(3 -len(excludedParams))
		signatureString = signatureString[:len(signatureString) - 1]
		signatureString = signatureString + ('->(n)')

		if(len(excludedParams) == 3):
			self.phi, self.phiError = self.phiSinge(E, r1, r2)
		else:
			self.phiPoints = np.vectorize(self.phiSinge, excluded = excludedParams, signature = signatureString)(E, r1, r2)
			self.phi, self.phiError = self.phiPoints.T
		
		return self.phi



class WKB1(WKB):

	def __init__(self, U, m, rmin = 10, rmax = 20, Emin = 1, Emax = 16, Ncrit = 10000):
		super().__init__(U, m, rmin, rmax, Emin, Emax, Ncrit)
		self.criticalPointsNum = 1

		self.getMaxPoint()

	def getMaxPoint(self, rGuess = None):

		pointWiseDerivatives = derivative(self.vectU, self.rLin, dx = 1E-6*((self.right - self.rleft)/self.Ncrit) )

		rmaxGuess = self.rLin[np.min(np.abs(pointWiseDerivatives)) == np.abs(pointWiseDerivatives)][0]

		if(rGuess):
			rmaxGuess = rGuess

		self.rmax = fmin(self.invU, rmaxGuess)[0] 
		self.Umax = self.U(self.rmax)

		return self.rmax

	def getTurningPoints(self, E, r1Guess = None, r2Guess = None):

		excludedParams = np.array([])

		if(not hasattr(E, '__len__')):
			excludedParams = np.append([0], excludedParams)

		if(not hasattr(r1Guess, '__len__')):
			excludedParams = np.append([1], excludedParams)

		if(not hasattr(r2Guess, '__len__')):
			excludedParams = np.append([2], excludedParams)

		signatureString = '(),'*(3 -len(excludedParams))
		signatureString = signatureString[:len(signatureString) - 1]
		signatureString = signatureString + ('->(n)')

		if(len(excludedParams) == 3):
			self.r1, self.r2 = self.getTurningPointsSingle(E, r1Guess, r2Guess)
			return np.array([self.r1, self.r2])
		else:
			self.turningPoints = np.vectorize(self.getTurningPointsSingle, excluded = excludedParams, signature = signatureString)(E, r1Guess, r2Guess)
			self.r1, self.r2 = self.turningPoints.T
			return self.turningPoints


	def getTurningPointsSingle(self, E, r1Guess = None, r2Guess = None):

		if(r1Guess):
			r1 = r1Guess
		else:
			r1 = self.rleft

		if(r2Guess):
			r2 = r2Guess 
		else:
			r2 = self.right


		rTurningLin = np.linspace(r1, r2, self.Ncrit)
		rLeftLin 	= rTurningLin[rTurningLin < self.rmax]	 
		rRightLin	= rTurningLin[rTurningLin >= self.rmax]

		diffF = lambda r: np.abs(self.U(r) - E)

		diffUE = np.abs(self.vectU(rTurningLin) - E)

		diffUELeft  = diffUE[rTurningLin < self.rmax]
		diffUERight = diffUE[rTurningLin >= self.rmax]

		rLeftGuess  = rLeftLin[np.min(diffUELeft) == diffUELeft]
		rRightGuess = rRightLin[np.min(diffUERight) == diffUERight]

		r1 = fsolve(diffF, rLeftGuess)[0]
		r2 = fsolve(diffF, rRightGuess)[0]

		return np.array([r1, r2])

	def phi(self, E):

		self.E = E
		self.getTurningPoints(self.E)

		return super().phi(E, self.r1, self.r2)

def phi(E, m, U, Emin = 0, Emax = 30, N = 100):

	Elin = np.linspace(Emin, Emax, N)

	K = - (2/hbar)

	f = lambda r: K * np.sqrt((2*m)*(U(r) - E))

	zerosPotential = lambda r: U(r) - E

	zerosArray = np.unique(fsolve(zerosPotential, Elin))

	zerosArray = np.unique(zerosArray[(zerosArray <= Emax) & (zerosArray >= Emin)])

	#print(zerosArray)

	#print(fsolve(zerosPotential, Emin))
	#print(fsolve(zerosPotential, Emax))

	plt.plot(Elin, np.vectorize(U)(Elin), color = 'blue')
	plt.axhline(E, color = 'red')
	plt.show()
	

	rmin, rmax = np.sort(zerosArray)[0], np.sort(zerosArray)[-1]

	#print(rmin, rmax)

	return quad(f, rmin, rmax)[0]





