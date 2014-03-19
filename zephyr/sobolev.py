#!/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python

'''
tool to calculate optical depths in region and 
access if the P_escape formalism has shortcomings
'''

import sv_sub as sv
import geo 
from constants import *




def tau_sobolev (nu, nl, gu, gl, f, freq, dvds):

	tau = nl - ( (gl / gu) * nu )

	tau *= PI_E2_OVER_M * ( f /  (freq * dvds) )

	return tau


def P_escape(tau):

	if tau < 1e-6:
		escape = 1.

	elif tau < 10.0:
		escape = (1. - np.exp (-tau)) / tau

	else:
		escape = 1. / tau

	return escape





'''
fname_to_open = sys.argv[1]
sv.input_params(fname = fname_to_open)



#define grid using input file parameters and user specified zscale, otherwise set to rstar
if nogrid_read:
	zscale=geo.rstar/10.0
	nz=30
xscale=geo.sv_rmin


xcen,zcen = geo.get_grid(zscale, xscale, nz, nx, geo.wind_rmax, 1)

xup, zup = geo.get_grid(zscale, xscale, nz, nx, geo.wind_rmax, 0)
'''



'''
define ds as a 1cm shell which we are considering
vmax is in km/sec and sets the velocity gradient in km/s/cm
'''


def get_escapes (vmax):
	'''
	define ds as a 1cm shell which we are considering
	vmax is in km/sec and sets the velocity gradient in km/s/cm
	'''
	vmin = 0.0 		# velocity at origin, 0

	dvmin = 1.0e-20


	thetas = np.arange(0.0, 0.5005*np.pi, 0.005*np.pi)
	theta_res = np.arange(len(thetas))



	dvds_by_angle = ( ( np.sin (thetas) + dvmin ) * vmax )
	dvds_ave = np.sum (dvds_by_angle) / len(dvds_by_angle)


	P_aves = []
	P_fudges = []

	Psum = 0
	tausum = 0



	for j in theta_res:

		dvds = dvds_by_angle[j]

		Psum += P_escape ( 1.0 / dvds )

		tausum += 1.0 / dvds 

	P_ave = Psum / len(theta_res)

	P_fudge = P_escape ( 1.0 / dvds_ave )

	tau_ave = tausum / len(theta_res)

	tau_fudge = 1.0 / dvds_ave

	return P_fudge, P_ave, dvds_by_angle, tau_ave, tau_fudge


v_array = np.arange( -1,2, 0.05)

fudges = []
aves = []
ratios = []
dvs = []

tau_ave = []
tau_dvds_ave = []

for logv in v_array:

	v = 10.0**logv

	P_fudge, P_ave, dvds_by_angle, tau, tauf= get_escapes (v)
	fudges.append(P_fudge)
	aves.append(P_ave)

	ratios.append(P_ave / P_fudge)

	#meandvds = np.sum(dvds_by_angle) / len(dvds_by_angle)

	tau_dvds_ave.append( tauf )
	tau_ave.append( tau )






import pylab


pylab.subplot(311)
pylab.plot( v_array, np.log10(fudges), label = "w/ dvds_ave")
pylab.plot( v_array, np.log10(aves), label = "P_ave" )
pylab.legend()

pylab.ylabel("P_escape")
pylab.ylim (-2, 0)


pylab.subplot(312)
pylab.plot(v_array, ratios)

pylab.xlabel("Log(dv/ds_max)")
pylab.ylabel("P_ave / P_dvdsave")



#pylab.subplot(313)
#pylab.subplot(313)
#pylab.plot( v_array, tau_ave , "k--", label = "tau_ave")
#pylab.plot( v_array, tau_dvds_ave, label = "w/ dvds_ave" )
#pylab.legend()
#pylab.xlabel("Log(dv/ds_max)")



pylab.show()





































