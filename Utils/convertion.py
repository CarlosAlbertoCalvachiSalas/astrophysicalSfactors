from sfactors.Utils.constants import *

def amuToMeVc_2(m):
	equivalent = physical_constants['atomic mass constant energy equivalent in MeV'][0]
	return m*equivalent

def amuToMeVfm_2s2(m):
	return (amuToMeVc_2(m)/(c**2))

