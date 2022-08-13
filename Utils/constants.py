from scipy.constants import *
import numpy as np

# SI UNITS 

k_ = 1.380649E-23
c_ = 299792458
alpha = 7.2973525693E-3
e_ 	= 1.602176634E-19
h_ = 6.62607015E-34 
hbar_ = h/(2*np.pi)

Dmass_ =	3.3435837724E-27
pmass_ =	1.67262192369E-27

# Nuclear units
# Energy MeV, distance fm, cross-secetions b

h = physical_constants['Planck constant in eV/Hz'][0]*1E-6
hbar = h/(2*np.pi)

c = c_ * 1E15

e = np.sqrt(alpha * (hbar * c))
