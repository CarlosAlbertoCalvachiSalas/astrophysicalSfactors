from sfactors.Solver.wkb import *
from sfactors.Solver.potentials import *
from sfactors.Databases.Reader.nuclearData import Nucleus

from sfactors.nonResonant.yakovlev import *
from testYakovlev import *

Elin = np.linspace(0.1, 16.26, 100)

wkbEx = WKB1(potential, Mg46Y.mu)
phiWKB = wkbEx.phi(Elin)
phiY = np.vectorize(Mg46Y.psi)(Elin) - (2*np.pi*np.vectorize(Mg46Y.eta)(Elin))

plt.plot(Elin, phiWKB, color = 'red', label = 'phi-WKB')
plt.plot(Elin, phiY, color = 'blue', label = 'phi-Y')

plt.legend()
plt.xlabel('E(MeV)')
plt.ylabel(r'$\Phi(E)$')
plt.title('WKB - exact yakovlev model calculation')
plt.show()

plt.plot(Elin, np.abs(phiWKB - phiY), color = 'black')
plt.xlabel('E(MeV)')
plt.ylabel(r'$|\Delta \Phi(E)|$')
plt.title('Residuals WKB - exact yakovlev model calculation')
plt.yscale('log')
plt.tight_layout()
plt.show()

rAtom = 1.16*((16)**(1/3))

wkbexpSingle = WKB1(lambda r: coulomb(r, 8, 8, rAtom) + expSingle(r, V0 = 1000, alpha = 0.5), Mg46Y.mu, rmin = 1)
wkbexpDouble = WKB1(lambda r: coulomb(r, 8, 8, rAtom ) + expDouble(r, V0 = 1/2000, alpha0 = 1, V1 = 1/2000, alpha1 = 1) , Mg46Y.mu, rmin = 1)
wkbWoodsSaxon = WKB1(lambda r: coulomb(r, 8, 8, rAtom ) + woodsSaxon(r, V0 = 1000, R = 1, a = 3), Mg46Y.mu, rmin = 1)
wkbGaussian = WKB1(lambda r: coulomb(r, 8, 8, rAtom ) + gaussian(r, V0 = 1000, alpha = 0.1), Mg46Y.mu, rmin = 1)

wkbY = WKB1(potential, Mg46Y.mu, rmin = 1)

r = np.linspace(0, 30, 1000)

plt.title(r'Coulomb potential charged sphere $V_c$')
plt.plot(r, np.vectorize(coulomb)(r, 8, 8), color = 'blue', label = 'coulomb-potential')
plt.ylabel('U(MeV)')
plt.xlabel('r(fm)')
plt.show()

wkbY.graphPotential(color = 'blue', label = 'yakovlev')
plt.title(r'Potential plot yakovlev')
plt.show()

wkbexpSingle.graphPotential(color = 'green', label = 'expSingle + coulomb')
plt.title(r'Potential plot $-V_{0}\exp{(-\alpha r)} + V_c $')
plt.show()

wkbexpDouble.graphPotential(color = 'orange', label = 'expDouble + coulomb')
plt.title(r'Potential plot $-\left(V_{0}\exp{(\alpha_0 r)} + V_{1}\exp{(\alpha_1 r)} \right)^{-1} + V_c$')
plt.show()

wkbWoodsSaxon.graphPotential(color = 'red', label = 'WS + coulomb')
plt.title(r'Potential plot $-V_{0}\left(1 + \exp{(r-R/a)} \right)^{-1} + V_c$')
plt.show()

wkbGaussian.graphPotential(color = 'black', label = 'gaussian + coulomb')
plt.title(r'Potential plot $-V_{0}\exp{(-\alpha r^2)} + V_c$')
plt.show()

#wkbWoodsSaxon = WKB1(lambda r: coulomb(r, 12, 12, rAtom ) + woodsSaxon(r, V0 = 1000, R = 1, a = 3), Mg46Y.mu, rmin = 1)
#phiWoodsSaxon = wkbWoodsSaxon.phi(Elin)
#phiY = wkbY.phi(Elin)

#plt.title('Comparison yakovlev-WS')
#plt.plot(Elin, phiY, color = 'blue', label = 'phi-Y')
#plt.plot(Elin, phiWoodsSaxon, color = 'red', label = 'phi-WS')
#plt.legend()

plt.show()

"""

wkbWoodsSaxon.graphPotential(color = 'red', label = 'WS')
wkbY.graphPotential(color = 'blue', label = 'Y')
wkbGaussian.graphPotential(color = 'black', label = 'Gaussian') 
wkbexpSingle.graphPotential(color = 'green', label = 'expSingle')
wkbexpDouble.graphPotential(color = 'orange', label = 'expDouble')

plt.show()

"""

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


