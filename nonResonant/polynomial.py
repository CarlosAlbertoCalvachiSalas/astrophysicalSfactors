import numpy as np
import pandas as pd 

def polynomial1(E, S0, g1):
	return S0 + g1*E

def polynomial2(E, S0, g1, g2):
	return S0 + g1*E + g2*(E**2)

def polynomial3(E, S0, g1, g2, g3):
	return S0 + g1*E + g2*(E**2) + g3*(E**3)

def polynomial4(E, S0, g1, g2, g3, g4):
	return S0 + g1*E + g2*(E**2) + g3*(E**3) + g3*(E**4)