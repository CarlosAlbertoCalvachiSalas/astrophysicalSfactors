from sfactors.nonResonant.yakovlev import *

import numpy as np
import matplotlib.pyplot as plt


alphaConstant = alpha(12, 12)
Ec = 16.27
delta =  0.0332
Rc = Rc(Ec = Ec, delta = delta, alpha = alphaConstant)
Rc1 = Rc1(Rc, delta)
beta = beta(delta)
Ec1 = Ec1(Ec, delta)

print(Rc, Rc1, beta, Ec1, sep = '\n')

RC0 = Rc*(1 - (beta**(-0.5)))
print(RC0)

def potential(r):

	#U = (alphaConstant/r)*(r >= Rc1) + (Ec * (1 - beta*(((r - Rc)**2)/(Rc**2))))*(r < Rc1)

	#return U * (U > 0)


	if(r >= Rc1):
		U = alphaConstant/r
	else:
		U = Ec * (1 - beta*(((r - Rc)**2)/(Rc**2)))

	if(U > 0):
		return U
	else:
		return 0.0


Mg46Y = Yakovlev(12, 12, 46, 46, 5.9785, 16.269, 0.0332, 0.5263 + (46 + 46)*(0.1393))


r = np.linspace(5, 20, 1000000)
U = np.vectorize(potential, otypes = [float])(r)
print(r[U > 0][U[U > 0] == np.min(U[U > 0])][0])

plt.plot(r, U, color = 'blue')
plt.xlabel('r(fm)')
plt.ylabel('U(MeV)')
#plt.show()

Mg46Y.plotSFactor()
#plt.show()

plt.close()
