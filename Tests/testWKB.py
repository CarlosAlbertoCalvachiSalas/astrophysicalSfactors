from sfactors.Solver.wkb import *
from sfactors.nonResonant.yakovlev import *
from testYakovlev import *


wkbEx = WKB1(potential, Mg46Y.mu)
wkbEx.graphPotential()

plt.show()

Elin = np.linspace(0.1, 16.26, 1000)

phiWKB = wkbEx.phi(Elin)
phiY = np.vectorize(Mg46Y.psi)(Elin) - (2*np.pi*np.vectorize(Mg46Y.eta)(Elin))


plt.plot(Elin, phiWKB, color = 'red', label = 'phi-WKB')
plt.plot(Elin, phiY, color = 'blue', label = 'phi-Y')

plt.legend()
plt.show()
#

#rleft, rright = wkbEx.getTurningPoints(E = Elin).T


#plt.plot(Elin, rleft, color = 'red')
#plt.show()


#plt.plot(Elin, rright, color = 'purple')
#plt.show()

#print(phi(10, Mg46Y.mu, np.vectorize(potential), Emin = 10, Emax = 25 ))

#phiVect = np.vectorize(phi, excluded = [1, 2])

#

#print('a')
#phiEval = phiVect(Elin, Mg46Y.mu, np.vectorize(potential))

#print(phiEval)

#plt.plot(Elin, phiEval, color = 'blue', label = 'Yakovlev-numerical')


