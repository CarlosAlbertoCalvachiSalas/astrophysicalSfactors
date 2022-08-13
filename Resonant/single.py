import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def breitWigner(E, Sr, Er, Gr):
	return Sr * (((Gr**2)/4))/((E-Er)**2 + ((Gr**2)/4))

def fitBreitWigner(E, S, p0 = None, sigma = None):
	return curve_fit(breitWigner, E, S, p0, sigma)

def generalResonance(E, Sr, Er, Gr, c):
	return Sr * (c)/((E-Er)**2 + ((Gr**2)/4))

