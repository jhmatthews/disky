
from astro_const import *

import numpy as np

def planck_nu (T, nu):

	''' 
	The planck function for frequency nu at temperature T
	'''
	x = H * nu / (BOLTZMANN * T)
	f = (2. * H * nu ** 3.) / (C ** 2. * (np.exp(x) - 1.))
    
	return f
	
	
	
	
def planck_lamb (T, lamb):

	''' 
	The planck function for wavelength lamb at temperature T
	'''
	x = H * C / (BOLTZMANN * T * lamb)
	f = (2. * H * C**2.0) / (lamb**5. * (np.exp(x) - 1.))
    
	return f
