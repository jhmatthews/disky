#! /Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
		recomb.py
		
Synopsis:
	An attempt to calculate level populations in an n level
	H atom, based only on recombination coefficients and 
	A values (Einstein A coefficient).

Usage:
	
Comments:
	This code uses
History:
	131023		Coding began (matom effort)
'''

from lines import A21, q12, q21, read_line_info
import numpy as np
from constants import *


def level_populations ( n, alphas, ne, line ):

	'''
	Calculates level populations for an n level atom populated only by recombination
	and cascades from upper levels.
    
    :INPUT:  
            n:  		int      
                		number of levels in atom
            alphas:	float array
            			recombination coefficients for levels.
            ne:		float
            			electron density

    :OUTPUT:
            levels:	array
            			level populations in number / cm**3
    :EXAMPLE:
            levels = level_populations ( n, alphas, ne )
            
    :COMMENTS:
    		careful with indices here. 0 in the alphas array is 1 in other arrays. 
    		
    		This is works out level populations by doing:
    			nrecombs direct to this level + number of transitions from upper levels to this level,
    			divided by the sum of the A coefficients for downwards transitions 
	'''
	
	npops = np.zeros(n)
	emiss = np.zeros(n)
	
	# first check you have recombination coefficients for levels in question
	if len(alphas) < n:
		print "Error: level_populations: Only %i alphas and a %i level atom" % (len(alphas),  n)
		return (-1)

	print "\n\t%s\t|\t%s\t|\t%s\t\t|\t%s\t" % ("Level", "Asum*n_i", "Asum", "nrecombs")
	print "------------------------------------------------------------------------------"
	
	# now cycle over levels, with the highest first	
	for i in range(n, 1, -1):
		
		# initially set n_i to be the number of recombinations direct to level i
		n_i = (ne * ne * alphas[i-1])

		# Calculating cascades into level- do all levels above i up to n
		# we start with the higher levels so as to get the correct populations
		# for lower levels populated from above
		for upper in range(i+1, n+1):
		
			for i_line in range(len(line)):
				
				if line[i_line].lu == upper and line[i_line].ll == i:

						#add the cascades from upper into i
						n_i += A21 ( line[i_line] ) * npops[upper-1]		# add the contribution
		
		Asum = 0
		



		#now sum up the A coefficients for all downward transitions from level i
		for lower in range(i-1, 0, -1):
		
			for i_line in range(len(line)):

				if line[i_line].lu == i and line[i_line].ll == lower:
				
					Asum += A21 ( line[i_line] )

		# dividing by the sum of A coefficients out of level i gives level populations
		n_i = n_i / Asum
		
		npops[i-1] = n_i
		
		
		# we know have level populations, can work out level emissivities = A_ij n_i h nu_ij
		for lower in range(i-1, 0, -1):
		
			emiss_sum = 0.
			for i_line in range(len(line)):

				if line[i_line].lu == i and line[i_line].ll == lower:
					emiss_sum += A21 ( line[i_line] ) * n_i * H * line[i_line].freq
			
			emiss[i-1] = emiss_sum
		
		print "\t%i\t|\t%.4f  \t|\t%8.2e\t|\t%.4f\t" %(i, n_i*Asum, Asum,  (ne * ne*alphas[i-1]) )
	
	return npops, emiss
		

print "\nCalculating strength of recombination lines for Case A.\n\n"
	
#Let's work with a 4 level Hydrogen atom, for now
filename = "data/atomic_macro/h4_lines.py"	# atomic data file

# this next line simply reads line info from the file and places in a line class instance
line_info = read_line_info(filename)



nlevels = 4			# 4 level macro atom
T = 10000.0			# 10000K
ne = 1.4e6			# electron density
nprots = 1.4e6		# H+ density
V = 1.3e23			# volume of cell




# recombination coefficients from Osterbrock. 
# the subscript gives the name of the subshell
# note that entry 0 is level 1
alpha_S = np.array([ 1.58e-13, 2.34e-14, 7.81e-15, 3.59e-15 ])
alpha_P = np.array([0.0, 5.35e-14, 2.04e-14, 9.66e-15])
alpha_D = np.array([0.0, 0.0, 1.73e-14, 1.08e-14])
alpha_F = np.array([0.0, 0.0, 0.0, 5.54e-15])

# sum the subshells to give recombiantion coefficient for level n
alpha_sum = alpha_S + alpha_P + alpha_D + alpha_F


	
# Get the level populations 
# function described above
level_pops, emiss = level_populations ( nlevels, alpha_sum, ne, line_info)


print "\n  Level\t|\t  Pops     \t|  Emissivity rel to level 4"
print "------------------------------------------------------"
for i in range(1,len(level_pops)):
	print "  %i   \t|\t %8.2e \t|\t%.2f" % (i+1, level_pops[i], emiss[i]/emiss[3])


# find where Halpha and Hbeta are in the line list
for i_line in range(len(line_info)):

	if line_info[i_line].lu == 3 and line_info[i_line].ll == 2: 	# then it is H alpha
		where_H_alpha = i_line
		nu_H_alpha = line_info[i_line].freq
		
	if line_info[i_line].lu == 4 and  line_info[i_line].ll == 2: 	# then it is H beta
		nu_H_beta = line_info[i_line].freq
		where_H_beta = i_line


# calculate line strengths for Halpha and Hbeta
H_alpha = H * nu_H_alpha * A21(line_info[where_H_alpha]) * level_pops[2]
H_beta = H * nu_H_beta * A21(line_info[where_H_beta]) * level_pops[3]



# print out the relative emissivity of Halpha to Hbeta
print "\nH_alpha / H_beta  %f\n" % (H_alpha / H_beta)


# print out all other line strengths
for i_line in range(len(line_info)):
		I = H *  A21(line_info[i_line]) * level_pops[line_info[i_line].lu - 1] *line_info[i_line].freq / H_beta
		print "Line strength %i => %i:  %.2f" % (line_info[i_line].lu, line_info[i_line].ll, I)
		
print"\n"




