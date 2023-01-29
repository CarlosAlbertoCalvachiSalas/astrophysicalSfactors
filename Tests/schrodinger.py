import numpy as np
import matplotlib.pyplot as plt
from sfactors.Solver.schrodinger import *

def potential(x):

	return 20*(x < 100)

shEq = SchrodingerEquation(m = amuToMeVfm_2s2(1), U = potential)

x, y = shEq.solve(E = 10).t, shEq.solve(E = 10).y

plt.plot(x, y[0] - np.exp(-shEq.k*x), color = 'blue')
#plt.plot(shEq.xLin, psi - , color = 'blue')
#plt.plot(shEq.xLin, potential(shEq.xLin), color = 'red')
#plt.plot(shEq.xLin, np.sinh(shEq.k*shEq.xLin), color = 'green', linestyle='dashed')
#plt.plot(x, np.exp(-shEq.k*x), color = 'purple', linestyle='dashed')
#plt.axvline(shEq.T, color = 'red')
plt.show()

