#! /Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
disk.py contains functions for things like 
disk luminosities, eddington fraction and alpha_ox
'''

from constants import *


def alp_ox (L_X, L_O):
	'''
	Calculates alpha ox for a given 2kev and 2500A luminosity	
	'''
	
	alpha = 0.3838 * np.log10( L_X / L_O )
	
	return alpha
	
	
def Ledd (m):
	'''
	calculates eddington luminosity for a solar mass m
	'''
	m *= MSOL
	
	consts = (4.0 * PI * G * C * MPROT ) / THOMPSON
	
	L = consts * m
	
	return L
	
	

def mdot_from_edd ( edd_frac, m ):
	''' 
	calculates an accretion rate from an eddington fraction
	mass m
	'''
	
	L = Ledd (m)		# eddington luminosity
	
	mdot = L / ( edd_frac * (C ** 2) )
	
	mdot = mdot * (3600.0 * 24.0 * 365.25) / MSOL
	
	return mdot


def L_two (L_X, alpha):
	'''
	L_two calculates the monochromatic X-ray luminosity at 2Kev
	
	Arguments:
		L_X 		2-10kev luminosity in ergs
		alpha	power law slope of spectrum
	'''
	
	f2 = 2000.0 / HEV	# freq at 2 kev
	f10 = 10000.0 / HEV 	# freq at 10 kev
	
	
	const_agn= L_X / ((( (f10**( alpha + 1.))) -  (f2**( alpha + 1.0))) /(alpha + 1.0))
	
	L = const_agn*f2**alpha

	return L
	
	
	
def L_2500 ( mdot, mbh ):
	'''
	L_two calculates the monochromatic X-ray luminosity at 2Kev
	
	Arguments:
		m		black hole mass
		mdot 	accretion rate
	'''
	f1 = 1.0e14; f2 = 1.0e18
	
	rmin = 0.5 * Schwarz ( mbh*MSOL )		# gravitational radius
	rmax = 1.0e17				# standard for models
	
	f, s = spec_disk (f1,f2,mbh,mdot,rmin,rmax)
	
	nu_2500 = C / (2500.0 * ANGSTROM)
	
	nu_ref = f[0]
	n=0
	
	while nu_ref < nu_2500:
		n += 1
		nu_ref = f[n]
	
	L = 0.5 * ( s[n] + s[n-1] )

	
	return L
	
	
def L_bol ( mdot, mbh ):
	'''Calculate Bolometric luminosity of a disk
	
	Arguments:
		m			mass of cental object in msol
		mdot			accretion rate in msol / yr
	'''
	
	rmin = 6.0 * 0.5 * Schwarz ( mbh*MSOL )		# 6 * gravitational radius
	rmax = 1.0e17				# standard for models
	
	f1 = 1.0e14; f2 = 1.0e18
	freq, spec = spec_disk (f1,f2,mbh,mdot,rmin,rmax)
	
	df = freq[1] - freq[0]
	sum_spec =  spec[0] * df
	
	for i in range(1, len(freq) - 1 ):
		print sum_spec
		df = freq[i] - freq[i-1]
		sum_spec +=  spec[0] * df
		
	sum_spec += spec[-1] * df
		
	return sum_spec



def Schwarz(m):
	'''calculate Schwarzschild radius for mass m''' 
	return 2.0 * G * m / (C**2)


	
def spec_disk ( f1, f2, m, mdot, rmin, rmax, nfreq = 1000, nrings = 100):
	'''
	spec_disk creates arrays of frequency and monchromatic luminosity for a disk
	
	Arguments:
		f1, f2 		frequency limits
		m			mass of cental object in msol
		rmin, rmax	minimum and maximum radius
		mdot			accretion rate in msol / yr
		nfreq		number of frequency points [optional]
		nrings 		number of disk annuli [optional]
		
	Also requires:
		constants.py, functions teff and tdisk
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
	reference temperature of a disk 
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




	
	
	































