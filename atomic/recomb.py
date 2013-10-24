#! /Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
		recomb.py
		
Synopsis:
	An attempt to calculate level populations in an n level
	H atom, based only on recombination coefficients and 
	A values (Einstein A coefficient).

Usage:
	
Comments:

History:
	131023		Coding began (matom effort)
'''

from lines import A21, q12, q21, read_line_info
import numpy as np
from constants import *


#Let's work with a 4 level Hydrogen atom, for now
filename = "data/atomic_macro/h4_lines.py"

T = 10000.0

emiss_filename = "emiss"


# load emissivities into array
emiss = np.loadtxt(emiss_filename, comments = '#', unpack=True, dtype = 'float')

emissivities = emiss[2]
nlevel = emiss[0]
nlevels = len(nlevel) - 1

	
print nlevel	

absorbed = emiss[1]


# at 10000K
# units are cm**3/s
alpha_S = np.array([ 1.58e-13, 2.34e-14, 7.81e-15, 3.59e-15 ])
alpha_P = np.array([0.0, 5.35e-14, 2.04e-14, 9.66e-15])
alpha_D = np.array([0.0, 0.0, 1.73e-14, 1.08e-14])
alpha_F = np.array([0.0, 0.0, 0.0, 5.54e-15])

alpha_sum = alpha_S + alpha_P + alpha_D + alpha_F


line = read_line_info(filename)
line_transitions = np.zeros ( (5, 5) )


recombs_to_levels = np.zeros(len(nlevel))

ne = 1.4e6
np = 1.4e6
V = 1.3e23

total_n = ne * V

for i in range(nlevels):
	nrecombs_to_level = total_n * alpha_sum[i]		# in units of number per second
	recombs_to_levels[i] = nrecombs_to_level
	print nrecombs_to_level
	
# recombs to levels now contains the level populations as a result of
# direct recombinations to that level, but we also have to get the 
# contributions from downwards cascades to that level from upper levels
# we now do the cascade downwards

for upper in range(4,0, -1):
	n_upper = recombs_to_levels[upper - 1]
	
	for lower in range(upper-1,0, -1):
	
		for i in range(len(line)):
			if line[i].lu == upper and line[i].ll == lower:
			
				print "\n-----------------------------"
				print "Calculating transitions from  %i to %i\n" % (upper, lower)
				print "Level pop %i:\t %8.4e" % (upper, n_upper)
				
				A_coeff = A21(line[i])
				print "A_coeff:\t %8.4e" % (A_coeff)
				
				line_transitions[upper][lower] = A_coeff * n_upper 	# line_transitions now contains number of photons for that level
				line_transitions[upper][lower] *= H * line[i].freq
				
				nlower_contribution = A_coeff * n_upper
				
				recombs_to_levels[lower - 1] += nlower_contribution
	

				
for upper in range(4,0, -1):
	for lower in range(upper-1,0, -1):
		print "Contribution from transition %i to %i is: %8.4e" % (upper, lower, line_transitions[upper][lower])


'''	
for i in range(len(line)):
	print line[i].lu,  line[i].ll
	
where_list_A = 0
where_list_B = 0

nlines = len(line)

for n in range(nlevels):
	i = n+1
	sum_out = 0
	print "\n-----------------------------"
	print "LEVEL %i" % i
	for l in range(nlines):
		if line[l].lu == i:
			term = recombs_to_levels[n] * A21(line[l])
			sum_out += term
			print "Line Transition %i to %i: %8.4e" % (line[l].lu, line[l].ll, term)
	print "Level %i: %8.4e" % (i , sum_out)'''

		
		






