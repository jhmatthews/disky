import numpy as n
import disky_const as cts
'''
Routines to Calculate the Continuum spectrum of a Hydrogen slab
given a certain temperature and column density. The source function 
is a blackbody of a given temperature. The radiative trasnfer only 
considers bound free and free-free absorptions.

Example:

import hyslab
wave=numpy.arange(3000,9000,0.1)	# Wavelength array in angstroms
flux=hyslab.hspec(wave,12000,21.5)	# Flux of spectra at T=12,000K and log(n)=21.5


Nov 2013 - JVHS	- First version. Adapted from STSCI synphot. Details on every routine.
'''
def bnu(wave,temp):
	'''
	BNU -- Planck function
	input:
		wave	- wavelength in angstroms
		temp	- kelvin temperature
	output:
		bnu		- blackbody intensity (erg/cm2/s/hz/ster)
	-
	nov 2013 J.V. Hernandez - Python
	jun 1989 Heith Horne @ stsci - clean up old subroutine
	jul 1989 Dave Bazell -- SPP version
	'''
	c1 = 1.43883e8			# hc/k 
	c2 = 1.95722e5			#

	x=wave*temp
	if x <= 0.0:			#Check for faulty data
		return(0.0)

	x = c1/x
	if x < 1.0e-4:
		factor = 2.0 / ( x * ( x + 2. ) )
	else:
		if x < 85.0:
			factor = 1. / ( n.exp( x ) - 1. )
		else:
			return(0.0)
	x=x*temp/c2
	return (x * x * x * factor)



def ophyd(wave,temp,z):
	'''
	ophyd -- hydrogenic cross section
	 bound-free plus free-free cross-section of a hydrogenic ion.
	 formulae are from seaton (stimulated emission is omitted.!!!)
	
	 input:
		wave	- wavelength (angstroms)
		temp    - kelvin temperature
		z	    - atomic number (1 = hi,2 = heii, etc.)
	 output:
		ophyd	- cross section per ion (cm^2/ion)
	--
	 1989 may kdh @ stsci .. modify comments
	 jul 1989 Dave Bazell -- SPP version
	'''
	sig0 = 7.907e-18
	k = 3
	# Check for bad initial values and return if there are any.
	if wave <= 0. or temp <= 0. or z <= 0:
		return(0.0)
	zz = z*z
	beta = 157890.*zz/temp
	x = wave*zz/911.267			# lambda*R , R=2 * pi^2 * m * e^4/ h^3 * c

	q = x**(-1.0/3.0)
	qq = q*q
	g1 = 0.3458*q
	g2 = 0.03307*qq
	
	# sum continua of first few levels explicitly
	n1 = int(n.sqrt(x)) + 1
	n2 = n1 + k - 1
	suma = 0.0
	
	for i in n.arange(n1,n2+1,1):
		nn = i*i
		xg = x/nn - 0.5				# 
		gbf = 1. + g1*xg - g2*(xg*xg + 1.25)  #  Eq. 8.9 Gray
		suma = suma + n.exp(beta/nn)*gbf/(nn*i)
	
	# use continuum approximation for remaining levels
	xg = x/beta
	gff = 1. + g1*(xg + 0.5) - 2.*g2*(xg*(xg + 0.5) + 0.75)
	a = n2 + 0.5
	suma = suma + (n.exp(beta/(a*a)) - 1. + gff)/(beta + beta)
	#print sig0*x*x*x*suma*n.exp(-beta)/zz
	return (sig0*x*x*x*suma/zz*n.exp(-beta))

def stimcor(wave,temp):
	'''
	STIMCOR -- Evaluates stimulated emission factor 1 - exp(-hnu/kt)
	 input:
		wave	- wavelength (angstroms)
		temp	- kelvin temperature
	 output:
		stimcor	- (1 - exp(-hnu/kt))
	--
	 may 1989 kdh @ stsci - revise comments
	 jul 1989 kdh @ stsci - use cgs to get physical constants
	 jul 1989 Dave Bazell -- SPP version
	'''
	
	factor=6.6262e-27*2.997925e10*1e8/1.3806e-16

	if wave <= 0. or temp <= 0.:
	   return(1.0)
	else:
	   return (1. - n.exp( -factor/(temp*wave) ) )


def hydnu(wave,temp,colden):
	'''
	HYDNU -- LTE hydrogen slab
	
	  emission from hydrogen.  bound-free and free-free 
	
	 Input:
		WAVE	= Wavelength (Angstroms)
		TEMP	= Kelvin Temperature
		COLDEN	= column density (cm^-2) or log(column density)
	 Output:
		HYDNU	= LTE source function (erg/s/cm^2/Hz/ster)
	--
	 Nov 2013 JVHS - Python transfer
	 Nov 1985 KDH @ STScI
	 Mar 1989 KDH @ STScI - use column density
	 Jul 1989 SPP version by Dave Bazell
	'''
	# Check for incorrect input values
	if wave <= 0. or temp <= 0.:
		return()

	# Calculate source function aka, blackbody
	value = bnu( wave, temp)
	if (value <= 0. ):
		return()

	# opacity
	z = 1.0 			# atomic Number Hydrogen
	opacity = ophyd( wave, temp, z ) * stimcor( wave, temp )
	if opacity <= 0.:
		return()

	# optical depth
	if ( colden > 0. and colden < 80. ):
		tau = opacity * (10**colden)
	else:
		tau = opacity * colden
	# radiative transfer
	if ( tau < 1.e-3 ):								# Optically thin
		
		factor = tau * ( 1. - tau*(0.5 - tau/6.) )
	else:
		factor = 1. - n.exp(-tau) 					# Optically thick
	return (value * factor)


def hspec(wave,temp,colden):
	'''
	HSPEC -- Calculate hydrogen spectra

	Input:
		WAVE	= Wavelength array (Angstroms)
		TEMP	= Kelvin Temperature
		COLDEN	= column density (cm^-2) or log(column density)
	 Output:
		hspec	= LTE source spectra for the wavelength array (erg/s/cm^2/Hz/ster)
	Nov 2013 JVHS - Python transfer
	'''
	flux=[]
	for i in wave:
		flux.append(hydnu(i,temp,colden)*2997924580000000000/(i**2))
	return(n.array(flux))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def opy(wave,temp,z):
	'''
	BB -- Calculate blackbody spectra
	'''
	flux=[]
	for i in wave:
		flux.append(ophyd(i,temp,z))
	return(n.array(flux))

def stim(wave,temp):
	'''
	BB -- Calculate blackbody spectra
	'''
	flux=[]
	for i in wave:
		flux.append(stimcor(i,temp))
	return(n.array(flux))


def bb_nu (wave,temp):

	''' 
	The planck function for frequency nu at temperature T
	'''
	nu=2.997e18/wave
	return ((1.474501e-47*nu**3)/(n.exp((4.799216e-11*nu)/(temp))-1.0))

def bb_lam (wave,temp):

	''' 
	The planck function for frequency nu at temperature T
	'''
	return ((1.1904397e27/wave**5)/(n.exp((1.438769e8)/(wave*temp))-1.0))

def bb(wave,temp,nu=True,old=False):
	'''
	BB -- Calculate blackbody spectra
	'''
	flux=[]
	if nu == True and old == False:
		for i in wave:
			flux.append(bb_nu(i,temp))
	if nu == True and old == True:
		for i in wave:
			flux.append(bnu(i,temp))
	if nu ==False:
		for i in wave:
			flux.append(bb_lam(i,temp))

	return(n.array(flux))
