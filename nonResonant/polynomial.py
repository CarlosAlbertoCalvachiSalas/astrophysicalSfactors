import numpy as np
import pandas as pd 
from sfactors.nonResonant.screening import *


def polynomial(n = 1, evaluate = True):
	
	intArray = np.arange(n) + 1
	intStrArray = intArray.astype(str)
		
	paramsArray = np.char.add('g', intStrArray)	
	expEStrArray = np.char.add('E**', intStrArray)

	termArray = np.char.add(paramsArray, '*(')
	termArray = np.char.add(termArray, expEStrArray)
	termArray = np.char.add(termArray, ')')

	termArray = np.append(['S0'], termArray)
	paramsArray = np.append(['E','S0'], paramsArray)

	termStr   = ' + '.join(termArray) 
	paramsStr = ', '.join(paramsArray)

	finalStr  = 'lambda ' + paramsStr + ': ' + termStr
	
	if(evaluate):
		return eval(finalStr)
	else:
		return paramsStr, termStr

#polynomial(n = 0)(1, 2)
#polynomial(n = 1)(1, 2, 3)
#polynomial(n = 4)(3, 4, 5, 6, 5, 8)
#polynomial(n = 10)(*np.arange(12))

def polynomial1(E, S0, g1):
	return S0 + g1*E

def polynomial2(E, S0, g1, g2):
	return S0 + g1*E + g2*(E**2)

def polynomial3(E, S0, g1, g2, g3):
	return S0 + g1*E + g2*(E**2) + g3*(E**3)

def polynomial4(E, S0, g1, g2, g3, g4):
	return S0 + g1*E + g2*(E**2) + g3*(E**3) + g3*(E**4)