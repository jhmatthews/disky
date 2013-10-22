import sys
import numpy as np
import pylab
import matplotlib.pyplot as plt
import scipy.integrate 
import scipy.optimize
from collections import namedtuple
import geo
import astro_help as ah
import disk_sub as disk

RADIAN=57.29598
C=2.997925e10
MSOL=1.979e33
G=6.670e-8
YR=3.1556925e7
EPSILON=1e-6
PI=3.1416
STEFAN_BOLTZMANN=5.669e-5


def tdisk (m, mdot, r):
	t = 3. * G / (8. * PI * STEFAN_BOLTZMANN) * m * mdot / (r * r * r)
  	t = pow (t, 0.25)
  	return (t)

def teff (t, x):
	q = (1.e0 -  (x ** -0.5e0)) / (x * x * x);
      	q = t * (q ** 0.25e0);
	return (q)


def spec_disk (f1,f2,m,mdot,rmin,rmax):
	tref=tdisk(m, mdot, rmin)
	nfreq=(f2/f1)*100
	freq=np.linspace(f1,f2,nfreq)
	spec=np.empty(nfreq)
	dfreq=freq[1]-freq[0]
	rtemp=np.logspace(np.log10(rmin),np.log10(rmax),num=100)
	rdisk=[]
	for j in range(len(rtemp)-1):
		rdisk.append((rtemp[j]+rtemp[j+1])/2.0)
		r=rdisk[j]/rmin
		area=PI*(rtemp[j+1]*rtemp[j+1]-rtemp[j]*rtemp[j])
		t=(disk.teff(tref,r))
		for i in range(len(freq)):
			spec[i]=spec[i]+(ah.planck_nu(t,freq[i])*area*PI*2)
	return (freq,spec)
			

def spec_disk1 (f1,f2,m,mdot,rmin,rmax):
	tref=tdisk(m, mdot, rmin)
	nfreq=1000
	freq=np.logspace(np.log10(f1),np.log10(f2),nfreq)
	spec=np.empty(nfreq)
	dfreq=freq[1]-freq[0]
	rtemp=np.logspace(np.log10(rmin),np.log10(rmax),num=100)
	rdisk=[]
	for j in range(len(rtemp)-1):
		rdisk.append((rtemp[j]+rtemp[j+1])/2.0)
		r=rdisk[j]/rmin
		area=PI*(rtemp[j+1]*rtemp[j+1]-rtemp[j]*rtemp[j])
		t=(disk.teff(tref,r))
		for i in range(len(freq)-1):
			spec[i]=spec[i]+(ah.planck_nu(t,(freq[i+1]+freq[i])/2.0)*area*PI*2*(freq[i+1]-freq[i]))
	return (freq,spec)


def lnu_disk (f,m,mdot,rmin,rmax):
	tref=tdisk(m, mdot, rmin)
	rtemp=np.logspace(np.log10(rmin),np.log10(rmax),num=100)
	rdisk=[]
	lnu=0.0
	for j in range(len(rtemp)-1):
		rdisk.append((rtemp[j]+rtemp[j+1])/2.0)
		r=rdisk[j]/rmin
		area=PI*(rtemp[j+1]*rtemp[j+1]-rtemp[j]*rtemp[j])
		t=(disk.teff(tref,r))
		lnu=lnu+(ah.planck_nu(t,f)*area*PI*2.0)
	return (lnu)


def llamb_disk (lamb,m,mdot,rmin,rmax):
	tref=tdisk(m, mdot, rmin)
	rtemp=np.logspace(np.log10(rmin),np.log10(rmax),num=100)
	rdisk=[]
	llamb=0.0
	for j in range(len(rtemp)-1):
		rdisk.append((rtemp[j]+rtemp[j+1])/2.0)
		r=rdisk[j]/rmin
		area=PI*(rtemp[j+1]*rtemp[j+1]-rtemp[j]*rtemp[j])
		t=(disk.teff(tref,r))
		llamb=llamb+(ah.planck_lamb(t,lamb)*area*PI*2.0)
	return (llamb)



def spec_disk2 (f1,f2,m,mdot,rmin,rmax):
	tref=tdisk(m, mdot, rmin)
	nfreq=10
	f1a=10**float(int(np.log10(f1)))
	f2a=10**float(int(np.log10(f2))+1)
	nrange=int(np.log10((f2a/f1a)))
	freq=[]
	dfreq=[]
	ftemp=f1a
	df=f1a/nfreq
	for i in range(nrange):
		for j in range(nfreq*9):
			ftemp=ftemp+df
			if ftemp > f2:
				break
			if ftemp >= f1:
				freq.append(ftemp)
		df=df*10.0


	#print freq[0],freq[len(freq)-1]
		

	spec=np.zeros(len(freq))
	rtemp=np.logspace(np.log10(rmin),np.log10(rmax),num=100)
	rdisk=[]
	for j in range(len(rtemp)-1):
		rdisk.append((rtemp[j]+rtemp[j+1])/2.0)
		r=rdisk[j]/rmin
		area=PI*(rtemp[j+1]*rtemp[j+1]-rtemp[j]*rtemp[j])
		t=(disk.teff(tref,r))
		for i in range(len(freq)-1):
			spec[i]=spec[i]+(ah.planck_nu(t,freq[i])*area*PI*2)
	return (freq,spec)

