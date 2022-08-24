import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from sfactors.nonResonant.polynomial import polynomial
from sfactors.nonResonant.screening import screeningFactor

def exponential(n = 1, evaluate = True):
	
	paramsStr, termStr = polynomial(n = n, evaluate = False)
	termStr = 'np.exp(' + termStr + ')'

	finalStr  = 'lambda ' + paramsStr + ': ' + termStr

	if(evaluate):
		return eval(finalStr)
	else:
		return paramsStr, termStr

#exponential(n = 0)(3, 4)
#exponential(n = 1)(3, 4, 5)
#exponential(n = 4)(*np.arange(4 + 2))
#exponential(n = 10)(*np.arange(10 + 2))

def exponential1(E, S0, g1):
	return S0 * np.exp(g1*E)

def exponential2(E, S0, g1, g2):
	return S0 * np.exp(g1*E + g2*(E**2))

def exponential3(E, S0, g1, g2, g3):
	return S0 * np.exp(g1*E + g2*(E**2) + g3*(E**3))

def exponential4(E, S0, g1, g2, g3, g4):
	return S0 * np.exp(g1*E + g2*(E**2) + g3*(E**3) + g3*(E**4))




def exponentialScreening(n = 1, evaluate = True):

	paramsStr, termsStr = exponential(n = n, evaluate = False)


def exponential1Screening(E, S0, g1, S, Ue, eta):

	return exponential(E, S0, g1) + screeningFactor(E, S, Ue, eta)
