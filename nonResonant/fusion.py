import numpy as np


def empiricalFusion34(E, b1, b2, b3, c1, c2, c3, c4, Ec, D):
	return np.exp(b1 + b2*E + b3*(E**2) + (c1 + c2*E + c3*(E**2) + c4*(E**3))/(1+np.exp((Ec - E)/D)))