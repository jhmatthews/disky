#! /Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
		recomb_sub.py
		
Synopsis:
	Subroutines for an attempt to calculate level populations in an n level
	H atom, based only on recombination coefficients and 
	A values (Einstein A coefficient).
	
	

Usage:
	
Comments:

History:
	131023		Coding began (matom effort)
'''
import numpy as np
from lines import A21, q12, q21, read_line_info, read_level_info
from constants import *



def levels_from_probs( n, eprbs, jprbs, eprbsnorm, jprbsnorm):
	'''
	Calculates level populations for an n level atom by considering the jumping
	probabilties '''
	return 0
	


def subshell_pops ( n, alphas, ne, level, rad ):

	'''
	Calculates level populations for an n level atom populated only by recombination
	and cascades from upper levels.
    
    :INPUT:  
            n:  		int      
                		number of levels in atom
            alphas:	float array
            			recombination coefficients for levels, split by subshell.
            			in form [alpha_s, alpha_p, etc...]
            ne:		float
            			electron density
            level	object array
					array of chianti_level class instances
			rad		object array
					array of chianti_rad class instances
	:OUTPUT:
            levels:	array
            			level populations in number / cm**3			
            		
     
     '''
	print 'SUBSHELL POPS:'
	
	emiss_principle = np.zeros(n)
	
	# first check you have recombination coefficients for levels in question
	if len(alphas[0]) < n:
		print "Error: level_populations: Only %i alphas and a %i level atom" % (len(alphas),  n)
		return (-1)
	
	# work out which levels we are using. n here is the maximum principal quantum number we are using
	levels_used=[]
	for i in range ( len (level)):
		if level[i].n < n: 
			levels_used.append ( level[i])
	
	
	# total levels is the sum of all subshells used and so is the length of the levels_used array
	total_levels = len(levels_used)
	
	npops = np.zeros(total_levels)
	emiss = np.zeros(total_levels)
	
	# now cycle over levels, with the highest first	
	# for level 1s the i index is 0
	# for the highest level e.g. 4f the i index is total_levels - 1.
	for i in range ( total_levels - 1, -1, -1):
	
		subshell = levels_used[i].notation[1]			# subshell string, e.g. s, p
		n_level = levels_used[i].n						# principal quantum number of level
		
		# alpha_index helps us choose the recombination coefficient
		# according to the subshell we are working with
		if subshell == 's': alpha_index = 0
		if subshell == 'p': alpha_index = 1
		if subshell == 'd': alpha_index = 2
		if subshell == 'f': alpha_index = 3
		
		
		# initially set n_i to be the number of recombinations direct to level i
		n_i = (ne * ne * alphas[alpha_index][n_level-1])
		
		
		# now we need to loop over all higher levels and work out their contribution to the 
		# level population, given by A12 * n2
		for j in range(i, total_levels):
		
			upper = levels_used[j].notation		# LS coupling notation for upper subshell
			lower = levels_used[i].notation		# LS coupling notation for lower subshell
			
			# loop over all lines in the wgfa file
			for i_line in range(len(rad)):
				
				if rad[i_line].note_up == upper and rad[i_line].note_low == lower:

						#add the cascades from upper into i
						n_i += rad[i_line].A * npops[j]		# add the contribution
		
		Asum = 0
	
		#now sum up the A coefficients for all downward transitions from level i
		for j in range(i-1, 0, -1):
		
			upper = levels_used[i].notation		# LS coupling notation for upper subshell
			lower = levels_used[j].notation		# LS coupling notation for lower subshell
		
			for i_line in range(len(rad)):

				if rad[i_line].note_up == upper and rad[i_line].note_low == lower:
				
					Asum += rad[i_line].A

		# dividing by the sum of A coefficients out of level i gives level populations
		n_i = n_i / Asum
		
		npops[i] = n_i
		
		
		# we know have level populations, can work out level emissivities = A_ij n_i h nu_ij
		for j in range(i, -1, -1):
		
			emiss_sum = 0.
			for i_line in range(len(rad)):

				if rad[i_line].note_up == upper and rad[i_line].note_low == lower:
					emiss_sum += rad[i_line].A * n_i * H * rad[i_line].freq
			
			emiss[i] = emiss_sum
			emiss_principle[ n_level] += emiss_sum
		
		#print "\t%i\t|\t%.4f  \t|\t%8.2e\t|\t%.4f\t" %(i, n_i*Asum, Asum,  (ne * ne*alphas[i-1]) )
		
	
	
	
	print emiss_principle
	return npops, emiss, emiss_principle
	
	
	


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
