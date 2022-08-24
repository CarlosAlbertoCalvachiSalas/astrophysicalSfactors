from sfactors.Utils.convertion import *

def alpha(Z1, Z2):
	return Z1*Z2*(e**2)

def Er(mu, alpha):
	return ((alpha**2)*mu)/(2*(hbar**2))

def eta(E, Er):
	return np.sqrt(Er/E)

def sommerfeld(E, Z1, Z2, mu, convert = True):
	if(convert):
		mu = amuToMeVfm_2s2(mu)
	else:
		mu = mu
		
	return eta(E, Er(mu, alpha(Z1, Z2)))

def reducedMass(m1, m2):
	return (m1 * m2)/(m1 + m2)

def sfactorFromCrossSection(crossSection, E, Z1, Z2, A1, A2):

	m1, m2 = amuToMeVfm_2s2(A1), amuToMeVfm_2s2(A2)	
	mu = reducedMass(m1, m2)

	etaValues = eta(E, Er(mu, alpha(Z1, Z2)))

	return crossSection*E*np.exp(2*np.pi*etaValues)
