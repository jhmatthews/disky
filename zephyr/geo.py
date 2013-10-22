########################################################
#
#				JM/NSH 10/03/13: GEO.PY
#
#	gives geo structure and defines the grid for zephyr
#
########################################################

import numpy as np
import math

mstar=0.0
rstar=0.0
disk_mdot=0.0
wind_rmax=0.0
wind_mdot=0.0
sv_rmin=0.0
sv_rmax=0.0
sv_thetamin=0.0
sv_thetamax=0.0                     
sv_lambda=0.0
sv_v_infinity=0.0
sv_r_scale=0.0
sv_alpha=0.0   
mdot_norm = 0.0
sv_v_zero = 6e5
sv_gamma = 1.0
disk_radmax=0.0
xlog_scale=1.0e7 #set some default values for log scales
zlog_scale=1.0e7
QSO=0

def get_grid(zscale, xscale, nz, nx, rmax, cen_or_not):

	'''defines your grid in x and z coordinates- being the cen/upperboundaries of cells
	depending on the variable cen_or_not. if not 0: upperbound, if 0: centre'''

	xmaxscaled= rmax / xscale
	zmaxscaled= rmax / zscale
	dlogx = (np.log10(xmaxscaled)) / (nx - 3)
	dlogz = (np.log10(zmaxscaled)) / (nz - 3)
	x,z = np.zeros(nx), np.zeros(nz)

	for ix in range(nx):
		if cen_or_not!=0: 
			x[ix]= xscale * (10.0** (dlogx * (ix)) )	
		else:
			if ix>0: x[ix]= xscale * (10.0** (dlogx * (ix-1)) )	 	

	for iz in range(nz):
		if cen_or_not!=0:
			z[iz]= zscale * (10.0** (dlogz * (iz)) )	
		else:
			if ix>0: zscale * (10.0** (dlogz * (iz-1)) )	

	return x,z






		
	
