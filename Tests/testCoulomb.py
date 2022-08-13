import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from sfactors.Solver.coulomb import *

t, y =  coulombDiff(1, 0).t, coulombDiff(1, 0).y

plt.plot(t, y[0], color = 'green')
plt.show()

rho = np.linspace(0.05, 20, 1000)
eta = 1
l = 0

F = coulombF(rho, eta, l)
G = coulombG(rho, eta, l)

angle = rho - (np.pi*(l/2)) - eta*np.log(2*rho) + np.angle(gamma(l + 1 + (1j*eta)))


plt.plot(rho, F, color = 'blue', label = r'$F_{1, 0}$')
plt.plot(rho, G, color = 'red', label = r'$G_{1, 0}$')
plt.plot(rho, np.cos(angle) , color = 'purple')
plt.plot(rho, np.sin(angle), color = 'green')
plt.legend()
plt.title('Repulsive')
plt.show()

eta = -1

F = coulombF(rho, eta, l)
G = coulombG(rho, eta, l)

plt.plot(rho, F, color = 'blue', label = r'$F_{1, 0}$')
plt.plot(rho, G, color = 'red', label = r'$G_{1, 0}$')
plt.legend()
plt.title('Attractive')
plt.show()