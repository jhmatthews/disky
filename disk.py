#! /Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
disk.py contains functions for things like 
disk luminosities, eddington fraction and alpha_ox
'''

from disky_const import *
import numpy as np


def alp_ox (L_X, L_O):
	'''
	Calculates alpha ox for a given 2kev and 2500A luminosity	
	'''
	
	alpha = 0.3838 * np.log10( L_X / L_O )
	
	return alpha
	
	
def Ledd (m):

	'''
	calculates eddington luminosity for a solar mass m
	
	Args:
		m	mass in solar masses
	
	returns:
		Eddington luminosity in ergs s^-1, float
	'''
	
	m *= MSOL
	
	consts = (4.0 * PI * G * C * MPROT ) / THOMPSON
	
	L = consts * m
	
	return L
	
	

def mdot_from_edd ( edd_frac, m ):

	''' 
	calculates an accretion rate from an eddington fraction
	Args:
		edd_frac		eddington fraction
		m			mass of central object in solar masses
	
	returns:
		mdot in solar masses / yr
	'''
	
	L = Ledd (m)		# eddington luminosity
	
	mdot = L / ( edd_frac * (C ** 2) )
	
	mdot = mdot * ( YEAR ) / MSOL	# normalise units
	
	return mdot


def L_two (L_X, alpha):

	'''
	L_two calculates the monochromatic X-ray luminosity at 2Kev
	
	Arguments:
		L_X 		2-10kev luminosity in ergs
		alpha	power law slope of spectrum
		
	Returns:
		monochromatic luminosity in units of erg /s /Hz
	'''
	
	f2 = 2000.0 / HEV	# freq at 2 kev
	f10 = 10000.0 / HEV 	# freq at 10 kev
	
	
	const_agn= L_X / ((( (f10**( alpha + 1.))) -  (f2**( alpha + 1.0))) /(alpha + 1.0))
	
	L = const_agn*f2**alpha

	return L
	
	
	
def L_2500 ( mdot, mbh, f1 = 1.0e14, f2 = 1.0e18	 ):

	'''
	L_2500 calculates the monochromatic luminosity at 2500 Angstroms
	
	Arguments:
		m		mass of cental object, solar masses
		mdot 	accretion rate, solar masses / yr
		f1, f2 	frequency bounds in Hz
	
	Returns:
		monochromatic luminosity in units of erg /s /Hz
	'''
		
	
	rmin = 6.0 * 0.5 * Schwarz ( mbh*MSOL )		# 6.0 * gravitational radius (ISCO)
	rmax = 1.0e17								# standard for JM models
	
	f, s = spec_disk ( f1, f2, mbh, mdot, rmin, rmax)
	
	nu_2500 = C / (2500.0 * ANGSTROM)
	
	nu_ref = f[0]
	n = 0
	
	while nu_ref < nu_2500:
		n += 1
		nu_ref = f[n]
	
	L = 0.5 * ( s[n] + s[n-1] )

	return L
	
	
def L_bol ( mdot, mbh ):

	'''
	L_bol calculates bolometric luminosity of a disk
	
	Arguments:
		m		mass of cental object, solar masses
		mdot 	accretion rate, solar masses / yr
		
	Returns:
		L_bol in units of erg /s
	'''
	
	rmin = 6.0 * 0.5 * Schwarz ( mbh )		# 6 * gravitational radius
	rmax = 1.0e17				# standard for models
	
	f1 = 1.0e14; f2 = 1.0e18
	freq, spec = spec_disk (f1,f2,mbh,mdot,rmin,rmax)
	
	# spec contains monochromatic luminosity do need to multiply by df
	df = freq[1] - freq[0]
	sum_spec =  spec[0] * df
	
	for i in range(1, len(freq) - 1 ):

		df = freq[i] - freq[i-1]
		sum_spec +=  spec[0] * df
		
	sum_spec += spec[-1] * df
		
	return sum_spec





def Schwarz (m):

	'''
	calculate Schwarzschild radius for mass m in solar masses.
	
	Arguments:
		m	mass in solar masses 
	Returns: 
		radius in cm
	''' 
	
	m *= MSOL
	return 2.0 * G * m / ( C**2 )






	
def spec_disk ( f1, f2, m, mdot, rmin, rmax, nfreq = 1000, nrings = 100):

	'''
	spec_disk creates arrays of frequency and monchromatic luminosity for a disk
	
	Arguments:
		f1, f2 		frequency limits Hz
		m			mass of cental object in msol
		rmin, rmax	minimum and maximum radius in cm
		mdot			accretion rate in msol / yr
		nfreq		number of frequency points [optional]
		nrings 		number of disk annuli [optional]
		
	Returns:
		spec, freq	2 arrays, one containing monchromatic luminosity (erg /s /cm**2 and 
					one containing freq in Hz
	'''
		

	# reference temperature of the disk
	mdot = mdot * MSOL; m = m * MSOL
	tref=tdisk(m, mdot, rmin)
	
	
	# number of frequencies specified as optional arguments, linear spaced array
	freq=np.linspace( f1, f2, nfreq)
	
	spec = np.empty(nfreq)
	dfreq = freq[1]-freq[0]
	
	# logarithmically spaced radii
	rtemp = np.logspace(np.log10(rmin), np.log10(rmax), num = nrings)
	
	rdisk = []
	
	# loop over annuli
	for j in range(len(rtemp)-1):
		
		# rdisk contains midpoint values for each annulus
		rdisk.append((rtemp[j]+rtemp[j+1])/2.0)
		
		# divide by min radius
		r =rdisk[j]/rmin
		
		# area of annulus
		area = PI * (rtemp[j+1]*rtemp[j+1] - rtemp[j]*rtemp[j])
		
		t = ( teff(tref,r) )		# effective temperature of annulus
		
		for i in range(len(freq)):

			spec[i] = spec[i] + ( planck_nu(t,freq[i]) * area * PI * 2.)
						
	return freq,spec
	




def tdisk (m, mdot, r):

	''' 
	tdisk gives the reference temperature of a disk 
	m	black holes mass
	r	minimum radius
	mdot accretion rate
	'''
	
	t = 3. * G / (8. * PI * STEFAN_BOLTZMANN) * m * mdot / (r * r * r)
  	t = pow (t, 0.25)
  	
  	return (t)

def teff (t, x):

	''' 
	effective temperature of a disk at a point x 
	t	reference temperature of the disk
	r	radius / minimum radius
	'''
	
	q = (1.e0 -  (x ** -0.5e0)) / (x * x * x)
	q = t * (q ** 0.25e0 )
	
	return (q)


def planck_nu (T, nu):

	''' 
	The planck function for frequency nu at temperature T
	'''
	x = H * nu / (BOLTZMANN * T)
	f = (2. * H * nu ** 3.) / (C ** 2. * (np.exp(x) - 1.))
    
	return f




	
	
	































