import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from sfactors import FILE_PATH
from sfactors.Utils.constants import * 
from sfactors.Utils.convertion import * 
from sfactors.Utils.formulas import * 
import os

params = pd.read_csv(FILE_PATH + 'Databases/ModelData/yakovlev.csv')

# Potential and astrophysical S-factor

def potential(r, alpha, Ec, Rc, Rc1):

	if(r >= Rc1):
		U = alpha/r
	else:
		U = Ec * (1 - beta*(((r - Rc)**2)/(Rc**2)))

	if(U > 0):
		return U
	else:
		return 0.0

def psi(E, Ec, Er, delta, beta, Ec1):

	xr = E/Ec1

	xl = delta*np.sqrt(beta*(Ec/(Ec - E)))
	
	gamma = (((2+3*delta)**(3/2))*np.sqrt(delta))/((1+delta)**2)

	if(E < Ec1):

		psir = 4*np.sqrt(Er/E)*(np.arcsin(np.sqrt(xr)) + np.sqrt(xr*(1-xr)) )

		psil = -gamma*np.sqrt(Er/Ec)*((Ec - E)/(Ec))*((np.pi/2) + np.arcsin(xl) + xl*np.sqrt(1 - (xl ** 2)))

		return psir + psil
	else:
		return np.pi*(2*np.sqrt(Er/E) - (gamma*((Ec - E)/(Ec))*np.sqrt(Er/Ec) ))

def sFactor(E, S0, Ec, xi, eta, psi):

	if(E < Ec):
		return S0 * np.exp(psi(E))
	else:
		return S0 * np.exp(2*np.pi*eta(E)) * np.sqrt(E/Ec) * (1 + xi*((E - Ec)/(E)))

def Rc(Ec, delta, alpha):
	return (alpha*(2+3*delta))/(2*Ec*(1+delta)**2)

def Rc1(Rc, delta):
	return Rc*(1+delta)

def beta(delta):
	return 1/(delta*(2+3*delta))

def Ec1(Ec, delta):
	return Ec * ((2+2*delta)/(2+3*delta))

class Yakovlev():

	def __init__(self, Z1, Z2, A1, A2, *args, **kwargs):
		self.alpha = alpha(Z1, Z2)
		self.mu    = amuToMeVfm_2s2(reducedMass(A1, A2))	
		self.Er    = Er(self.mu, self.alpha)
		
		self.database = params[(params['Z1'] == Z1) & (params['Z2'] == Z2)].iloc[0]
		
		self.eta   = lambda E: eta(E, self.Er)
		self.S0 	= self.database['S0']

		dr1 = self.database['R1a'] if A1 >= 2*Z1 else self.database['R1b']
		dr2 = self.database['R2a'] if A2 >= 2*Z2 else self.database['R2b']
		rc0 = self.database['R'] + dr1*(A1 - 2*Z1) + dr2*(A2 - 2*Z2)	 
		
		self.Ec 	= self.alpha/rc0
		self.delta 	= self.database['delta']
		self.xi 	= self.database['xi0'] + self.database['xi1']*(A1 + A2)

		if(type(self.Ec) == type(None)):
			self.Ec = self.alpha/(R + R1*np.abs(A1 - 2*Z1) + R2*np.abs(A2 - 2*Z2))
			print('Ec = ', self.Ec, ' MeV')

		self.Rc   = Rc(self.Ec, self.delta, self.alpha)
		self.Rc1  = Rc1(self.Rc, self.delta)
		self.beta = beta(self.delta)
		self.Ec1  = Ec1(self.Ec, self.delta)

		self.potential = lambda r: potential(r, self.alpha, self.Ec, self.Rc, self.Rc1)

		self.psi = lambda E: psi(E, self.Ec, self.Er, self.delta, self.beta, self.Ec1)

		self.sFactor = lambda E: sFactor(E, self.S0, self.Ec, self.xi, self.eta, self.psi)

	def plotSFactor(self):

		Elin = np.linspace(0.005, 30, 100)
		sFactorArray = np.vectorize(self.sFactor)(Elin)

		plt.plot(Elin, sFactorArray, color = 'blue')

		plt.yscale('log')

		#plt.show()
# Useful formulas


def toEV(E):
	return E/e

def toJFromEV(E):
	return E*e


# Ueda functions

def sommerfeldParam(E, Z1, Z2, mu):
	return alpha*Z1*Z2*np.sqrt((mu*(c**2))/(2*E))

def gamowEnergy(Z1, Z2, mu):
	return 2*mu*(c**2)*((2*np.pi*Z1*Z2*alpha)**2)

# Fit functions



# From Ueda Paper

"""
def polynomial(dk):
	if(dk == 0):
		return lambda n: 1
	elif(dk == 2):
		return lambda n: 12*(n**2) + 18*n + 5
	elif(dk == 4):
		return lambda n: 144*(n**4) + 336*(n**3) + 81*(n**2) - 144*n -35
	elif(dk == 6):
		return lambda n:  1728*(n**6) + 4320(n**5) - 4320(n**4) - 13320(n**3) - 288(n**2) + 6210*n + 665

"""
