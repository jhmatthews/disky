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
		
	# line is a class instance like the line_ptr in PYTHONRT
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
		
	






















	
