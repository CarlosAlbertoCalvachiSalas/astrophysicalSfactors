# bpm stands for Barrier Penetration Model
import numpy as np 

from sfactors.Utils.constants import *
from sfactors.Solver.potentials import coulomb
from sfactors.Solver.wkb import *

def centrifugal(r, l, mu):
	return (l*(l+1)*(hbar**2))/(2*mu*(r**2))

def Tl(E, l, U, mu):

	effPot = lambda r: centrifugal(r, l, mu) + U(r)

	wkbPhi = WKB1(effPot, mu).phi(E)

	return 1/(1 + np.exp(wkbPhi))

def crossSection(E, mu, U, lmax = 3):

	k = np.sqrt((2*E*mu)/(hbar**2))

	C = np.pi/(k**2)

	l = np.arange(lmax + 1)

	TlArray = np.vectorize(Tl)(E, l, U, mu)

	sumResult = np.sum((2*l + 1)*TlArray)

	return sumResult * C
	# TO be programmed