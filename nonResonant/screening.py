import numpy as np
import matplotlib.pyplot as plt
from sfactors.Utils.formulas import *
from sfactors.Databases.Reader.nuclearData import Nucleus

def screeningFactor(E, S, Ue, Z1, Z2, mu):
	
	eta = sommerfeld(E, Z1, Z2, mu)

	return S*(E/(E+Ue))*np.exp((np.pi*eta*Ue)/E)

def screenedFunction(E, f, S, Ue, Z1, Z2, mu):
	screening = screeningFactor(E, S, Ue, Z1, Z2, mu)

	return f(E)*screening

d = Nucleus(1, 2)
mu = d.m/2

Elin = np.linspace(2E-3, 40E-3, 10000)
Ue = 309E-6


plt.plot(Elin*1E3, screeningFactor(Elin, S = 40, Ue = Ue, Z1 = d.Z, Z2 = d.Z, mu = mu), color = 'blue')
plt.xscale('log')
plt.xlim((1, 40))
#plt.plot(Elin, eta(Elin), color = 'red')
plt.show()