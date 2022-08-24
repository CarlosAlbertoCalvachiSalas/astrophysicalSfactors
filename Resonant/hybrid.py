import numpy as np
from .single import *

from sfactors.nonResonant import *

def breitWignerBackground1(E, Sr, Er, Gr, S0, g1):
	return breitWigner(E, Sr, Er, Gr) + exponential1(E, S0, g1)

def rmatrixPolynomial1(E, Sr, Er, Gr, c, S0, g1):
	return generalResonance(E, Sr, Er, Gr, c) + polynomial1(E, S0, g1)

def breitWignerBackground2(E, Sr, Er, Gr, S0, g1, g2):
	return breitWigner(E, Sr, Er, Gr) + exponential2(E, S0, g1, g2)

def rmatrixPolynomial2(E, Sr, Er, Gr, c, S0, g1, g2):
	return generalResonance(E, Sr, Er, Gr, c) + polynomial2(E, S0, g1, g2)
