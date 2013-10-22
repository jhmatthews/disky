#!/usr/bin/env python 
########################################################
#
#				SV_SUB.PY
#
#	NSH subroutines related to Shlosman Vitello wind
#
########################################################

import sys
import numpy as np
import pylab
import matplotlib.pyplot as plt
import scipy.integrate 
import scipy.optimize
from collections import namedtuple
import geo

RADIAN=57.29598
C=2.997925e10
MSOL=1.979e33
G=6.670e-8
YR=3.1556925e7
EPSILON=1e-6
PI=3.1416

def input_params(fname="sv68.pf"):

	geo.fname=fname
	inp=open(fname,'r')
	for line in inp.readlines():
		data=line.split()
		if (data[0][0:5]=="mstar") :
			geo.mstar=float(data[1])*MSOL 
		elif (data[0][0:5]=="rstar"):
			geo.rstar=float(data[1])
		elif (data[0][0:9]=="disk.mdot"):
			geo.disk_mdot=float(data[1])*MSOL/YR
		elif (data[0][0:16]=="QSO_BH_radiation"):
			geo.QSO=int(data[1])
		elif (data[0][0:7]=="lum_agn"):
			geo.lum_agn=float(data[1])
		elif (data[0][0:19]=="agn_power_law_index"):  
			geo.alpha_agn=float(data[1])
		elif (data[0][0:11]=="wind.radmax"):
			geo.wind_rmax=float(data[1])
		elif (data[0][0:9]=="wind.mdot"):
			geo.wind_mdot=float(data[1])*MSOL/YR
		elif (data[0][0:11]=="disk.radmax"):
			geo.disk_radmax=float(data[1])  
		elif (data[0][0:10]=="sv.diskmin"):
			geo.sv_rmin=float(data[1])*geo.rstar
		elif (data[0][0:10]=="sv.diskmax"):
			geo.sv_rmax=float(data[1])*geo.rstar
		elif (data[0][0:11]=="sv.thetamin"):
			geo.sv_thetamin=float(data[1])/RADIAN
		elif (data[0][0:11]=="sv.thetamax"):
			geo.sv_thetamax=float(data[1])/RADIAN                              
		elif (data[0][0:18]=="sv.mdot_r_exponent"):
			geo.sv_lambda=float(data[1])
		elif (data[0][0:13]=="sv.v_infinity"):
			geo.sv_v_infinity=float(data[1])
		elif (data[0][0:22]=="sv.acceleration_length"):
			geo.sv_r_scale=float(data[1])
		elif (data[0][0:24]=="sv.acceleration_exponent"):
			geo.sv_alpha=float(data[1])  
                      
	geo.sv_gamma=1.0
  	geo.mdot_norm =scipy.integrate.romberg (sv_wind_mdot_integral, geo.sv_rmin, geo.sv_rmax, rtol=1e-6);
	return (0)



def sv_velocity (x):
	v=np.zeros(3)
	zzz = v_escape = -99.
	rzero = sv_find_wind_rzero (x)
  	theta = sv_theta_wind (rzero)
 	r = np.sqrt(x[0] * x[0] + x[1] * x[1])
  	ldist = np.sqrt ((r - rzero) * (r - rzero) + x[2] * x[2])
 	vl = geo.sv_v_zero
  	if (ldist > 0):
    		zzz = pow (ldist / geo.sv_r_scale, geo.sv_alpha)
		if (rzero < geo.rstar):
			v_escape = np.sqrt (2. * G * geo.mstar / geo.rstar);
      		else:
			v_escape = np.sqrt (2. * G * geo.mstar / rzero);
		vl =geo.sv_v_zero + (geo.sv_v_infinity * v_escape -geo.sv_v_zero) * zzz / (1. + zzz)
	
	v[0] = vl * np.sin (theta)
	if (r > 0):
    		v[1] = np.sqrt (G * geo.mstar * rzero) / r
  	else:
    		v[1] = 0

  	v[2] = vl * np.cos (theta)

  	if (x[2] < 0):		
    		v[2] *= (-1)
  	if (x[1] != 0.0):
      		project_from_cyl_xyz (x, v, xtest)
     	 	v=stuff_v (xtest)
  	speed = (np.sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
  	return (speed,v)



def sv_rho (x):

	speed,v=sv_velocity (x)
	rzero = sv_find_wind_rzero (x)
  	theta = sv_theta_wind (rzero)
	r = np.sqrt (x[0] * x[0] + x[1] * x[1])
  	ldist = np.sqrt ((r - rzero) * (r - rzero) + x[2] * x[2])
	dmdot_da =geo.wind_mdot * (rzero**(geo.sv_lambda)) * np.cos (theta) / geo.mdot_norm / 2.
	if (rzero > geo.sv_rmax):
    		rzero = geo.sv_rmax;
  	dtheta_drzero =(sv_theta_wind (rzero)-sv_theta_wind ((1. - EPSILON) * rzero)) / (EPSILON * rzero)
	dr_drzero = 1. + ldist * dtheta_drzero / np.cos (theta)

  	rho = rzero * dmdot_da / (dr_drzero * r * v[2])
	return (rho)



def sv_find_wind_rzero (p):

  	z = np.fabs (p[2])
	if (z == 0):
    		x = (np.sqrt (p[0] * p[0] + p[1] * p[1]))
     		return (x)
	zero_p=sv_zero_init (p)	
  	rho_min = geo.sv_rmin + z * np.tan (geo.sv_thetamin)
  	rho_max = geo.sv_rmax + z * np.tan (geo.sv_thetamax)
  	rho = np.sqrt (p[0] * p[0] + p[1] * p[1])
  	if (rho <= rho_min):
    		x = geo.sv_rmin * rho / rho_min
      		return (x)
     	if (rho >= rho_max):
     		x = geo.sv_rmax + rho - rho_max
      		return (x)
    	x = scipy.optimize.brentq (sv_zero_r, geo.sv_rmin, geo.sv_rmax, args=zero_p,xtol=100.)
  	return (x)



def sv_zero_init (p):
  	zero_p= stuff_v (p)
  	zero_p[2] = np.fabs (zero_p[2])
  	return (zero_p)


def stuff_v (vin):
	vout=np.zeros(3)
  	vout[0] = vin[0];
  	vout[1] = vin[1];
  	vout[2] = vin[2];
  	return (vout)

def sv_theta_wind (r):
  	if (r <= geo.sv_rmin):
    		return (np.arctan (np.tan (geo.sv_thetamin * r / geo.sv_rmin)));
  	if (r >= geo.sv_rmax):
    		return (geo.sv_thetamax);
  	theta = geo.sv_thetamin + (geo.sv_thetamax - geo.sv_thetamin) *  ((r - geo.sv_rmin) / (geo.sv_rmax - geo.sv_rmin)**geo.sv_gamma)
  	return (theta);


def sv_wind_mdot_integral (r):
  	x = 2 * PI *  (r**( geo.sv_lambda + 1.)) * np.cos (sv_theta_wind (r))
  	return (x)

def sv_zero_r (r,zero_p):  
	theta = sv_theta_wind (r);
 	rho = np.sqrt (zero_p[0] * zero_p[0] + zero_p[1] * zero_p[1])
  	rho_guess = r + np.tan (theta) * zero_p[2]
  	return (rho_guess - rho)


