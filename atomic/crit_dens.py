

from lines import *


print '''Hi there -- This program calculates critical densities in Hydrogen. 
The critical density is when collisional dexcitation is negligible, 
and is given by the ratio of the sum over all lower levels of the 
Einstein A21 coefficient to the sum over all other levels of the 
collisional coefficients, q21 and q12.\n\n'''

filename = "data/atomic_macro/h20_lines.py"
line = read_line_info (filename)

nentries = len(line)


reference_array = np.arange(1,21)	# levels to look at.

n = len(reference_array)

'''
      omega =
	ECS_CONSTANT * line_ptr->gl * gaunt * line_ptr->f / line_ptr->freq;
      q21_a = 8.629e-6 / (sqrt (t) * line_ptr->gu) * omega;
      q21_t_old = t;'''



print "A21 for Halpha %8.4e" % A21(line[37])
print "q21 for Halpha %8.4e" % q21(line[37], 10000)
print "q12 for Halpha %8.4e" % q12(line[37], 10000)
print "g ratios for Halpha %.2f\n\n" % (1.0*line[19].g_u / line[19].g_l)


print " Level   |      Critical Density    |"
print "---------------------------------------"

critical_density(reference_array, 10000)

# All done
