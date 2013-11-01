#! /Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
	lines.py

Synopsis:
	This program calculates critical densities in Hydrogen.
	The critical density is when collisional dexcitation is negligible
	and is given by the ratio of the sum over all lower levels of the 
	Einstein A21 coefficient to the sum over all other levels of the 
	collisional coefficients, q21 and q12.

	This routine also contains functions q21, q12 and A21 for use elsewhere.

History:
	131017	Coding began

'''


# import constants and modules
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from constants import *
import classes as cls




def q21(line, T):

	'''calculate q21 for a line at temperature T'''
	omega = ECS_CONSTANT * line.g_l * line.osc / line.freq
	term = 8.629e-6 / (np.sqrt (T) * line.g_u) * omega;

	#print np.exp ( -H_OVER_K * line.freq / T ), np.sqrt(T)
	#print np.exp ( -H_OVER_K * line.freq / 20000.0 ), np.sqrt(20000.0)


	return term


def q12(line, T):

	'''calculate q12 for a line at temperature T'''

	term = (line.g_u / line.g_l) * q21 ( line, T) * np.exp ( -H_OVER_K * line.freq / T )

	return term



def A21(line):
	
	'''calculate A21 for a given line'''

	term_a = (line.freq**2) * line.osc
	term_b = (1.0*line.g_l) / (1.0*line.g_u)
	term_c = A21_CONSTANT

	return term_a * term_b * term_c 
	
	
	






def read_line_info(filename):

	'''read in line info from Python's atomic data e.g. h20_lines'''

	line_array_read = np.loadtxt(filename, comments='#', unpack = True, dtype = 'string')
	
	line_array_read = np.transpose(line_array_read)

	line = np.ndarray( len(line_array_read),dtype=np.object)
	
	for i in range(len(line)):
		z = float (line_array_read[i][1])
		ion = float (line_array_read[i][2])
		wave = ANGSTROM * float (line_array_read[i][3])
		freq = C / ( wave ) 
		osc = float (line_array_read[i][4])
		gl = int (line_array_read[i][5])
		gu = int (line_array_read[i][6])
		ll =  int (line_array_read[i][9])
		lu =  int (line_array_read[i][10])
		
		line[i] = cls.line(z, ion, wave, freq, osc, gl, gu, ll, lu)
		
	# line is an array of class instances like the line_ptr in PYTHONRT
	return line





def critical_density(reference_array, T):

	n = len(reference_array)
	filename = "data/atomic_macro/h20_lines.py"
	line = read_line_info (filename)
	for i in range(n):	# cycle over number of levels we want to look at 

		i_lev = reference_array[i]	# level to calculaye

		qsum = 0	# sum of q coefficients
		Asum = 0	# sum of A coefficients

		for j in range( len(line) ):	# go through line list

			if line[j].lu == i_lev:	# then we need to calculate A21 and q21

				#This means we have j < i, and thus need to calculate  
				#print 'j < i: ', line[j].lu, line[j].ll, i_lev, q21 ( line[j], T), A21 ( line[j])

				qsum += q21 ( line[j], T)
				Asum += A21 ( line[j])

			if line[j].ll == i_lev: 	# then we need to calculate q12, collisional excitation rate

				#print 'j != i: ', line[j].ll, line[j].lu, i_lev

				term = q12 ( line[j], T)
				#print term
				qsum += q12 ( line[j], T)


		crit_density = Asum / qsum	# critical density is the ratio of the 2 sums
		print "  %2i     |          %8.4e      |" % (i_lev, crit_density)
	return crit_density
		
	




def read_level_info(filename):

	'''read in level info from Python's atomic data e.g. h20_level'''

	level_array_read = np.loadtxt(filename, comments='#', unpack = True, dtype = 'string')
	
	level_array_read = np.transpose(level_array_read)

	level = np.ndarray( len(level_array_read),dtype=np.object)
	
	for i in range(len(level)):
		z = int (level_array_read[i][1])
		ion = int (level_array_read[i][2])
		lvl = int (level_array_read[i][3])
		ionpot = float (level_array_read[i][4]) 
		E = float (level_array_read[i][5])
		g = int (level_array_read[i][6])
		rad_rate = float (level_array_read[i][7])
		brack =  str (level_array_read[i][8])
		nnstring =  str  (level_array_read[i][9])
		
		level[i] = cls.level(z, ion, lvl, ionpot, E, g, rad_rate, brack, nnstring)
		
	# level is an array of class instances like the line_ptr in PYTHONRT
	return level


def sobolev(line, d1, d2, dvds):

	'''returns the sobolev optical depth for a line class instance
	
	:INPUT
		line		class instance
				information for the line from atomic data
		d1		float
				number density of lower level
		d2		float
				number density of lower level
		dvds		float
				velocity gradient
				
	:OUTPUT
		tau		float
				sobolev optical depth
			
	'''
	
	tau = (d1 - ((line.gl / line_ptr.gu) * d2 ) );
	tau *= PI_E2_OVER_M * line.f / line.freq / dvds;
	
	return tau




def read_chianti_data ( level_filename="h_1.clvlc", radiative_filename="h_1.wgfa"):

	'''
	read in a chianti atomic database and place in chianti class instance
	Default filenames are for hydrogen
	
	:INPUT
		level_filename		string
							filename for level information
		radiative_filename	string
							filename for radiative information
							
	:OUTPUT
		level	object array
				array of chianti_level class instances
		rad		object array
				array of chianti_rad class instances
	'''
	
	# read in data, place in array
	level_array_read = np.loadtxt (level_filename, comments = "%", dtype = "string")
	rad_array_read = np.loadtxt (radiative_filename, comments = "%", dtype = "string")
	
	# create blank arrays of type object to store class instances
	level = np.ndarray( len(level_array_read),dtype=np.object)
	rad = np.ndarray( len(rad_array_read),dtype=np.object)
	
	# create level class instances
	for i in range(len(level_array_read)):
		index = int (level_array_read[i][0])
		config = int (level_array_read[i][1])
		notation = str (level_array_read[i][2])
		spin = int (level_array_read[i][3]) 
		l = int (level_array_read[i][4])
		l_symbol = str (level_array_read[i][5])
		j = float (level_array_read[i][6])
		multiplicity = int (level_array_read[i][7])
		E_obs = float (level_array_read[i][8])
		E_obs2 = float (level_array_read[i][9])
		E_th = float (level_array_read[i][10])
		E_th2 = float (level_array_read[i][11])
		n = int(notation[0])
		
		level[i] = cls.level (index,  config, notation, spin, l, l_symbol, j, multiplicity, E_obs, E_obs2, E_th, E_th2, n)
	
	# create wgfa class instances
	for i in range(len(rad_array_read)):
		self.ll = int (rad_array_read[i][0])
		self.lu = int (rad_array_read[i][1])
		self.wave = float (rad_array_read[i][1])
		self.osc = float (rad_array_read[i][1])
		self.A = float (rad_array_read[i][1])
		
		rad[i] = cls.chianti_rad(ll, lu, wave, osc, A)
		
	return level, rad





	
