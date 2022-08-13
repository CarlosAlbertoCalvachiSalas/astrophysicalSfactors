import numpy as np
from .single import *

from sfactors.nonResonant import *

def breitWignerBackground(E, Sr, Er, Gr, S0, g1):
	return breitWigner(E, Sr, Er, Gr) + exponential1(E, S0, g1)

def rmatrixPolynomial(E, Sr, Er, Gr, c, S0, g1):
	return generalResonance(E, Sr, Er, Gr, c) + polynomial1(E, S0, g1)

