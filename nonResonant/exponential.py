import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

def exponential(order = 1):
	
	paramsArray = np.arange(order) + 1

	if(len(paramsArray) == 0):
		return lambda E, S0: S0
	else:

		strParamsArray = paramsArray.astype(np.str)
		paramsStr = ', '.join(np.char.add('g', strParamsArray))
		paramsStr = 'S0, ' + paramsStr

		print(paramsStr)

def exponential1(E, S0, g1):
	return S0 * np.exp(g1*E)

def exponential2(E, S0, g1, g2):
	return S0 * np.exp(g1*E + g2*(E**2))

def exponential3(E, S0, g1, g2, g3):
	return S0 * np.exp(g1*E + g2*(E**2) + g3*(E**3))

def exponential4(E, S0, g1, g2, g3, g4):
	return S0 * np.exp(g1*E + g2*(E**2) + g3*(E**3) + g3*(E**4))