import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 

from sfactors.Databases.Reader.nuclearData import *
from sfactors.Utils.convertion import amuToMeVc_2
from scipy.constants import *

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "sans-serif"
})

mp = physical_constants["proton mass in u"][0]
mn = physical_constants["neutron mass in u"][0]

B = mp*atomicNumbers + mn*(massNumbers - atomicNumbers) - atomicMass

Ba = amuToMeVc_2(B / massNumbers)

plt.scatter(massNumbers, Ba, s = 10, color = 'black')
plt.title('Binding energy per nucleon vs mass number curve')
plt.ylim((0, 9))
plt.ylabel(r'$B/A$ (MeV)')
plt.xlabel(r'$A$')
plt.savefig('BindingEnergyCurve.eps', format = 'eps')
plt.show()


