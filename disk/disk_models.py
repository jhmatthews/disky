from disky_const import *
from scipy import integrate, spatial
import numpy as np


class model:

	'''
	A class to store a set of stellar atmosphere models (or similar).
	To specify a stellar atmosphere model, all you need are T and g

	To initialise:

			mymodel = model("modelname")
			model.read()

		where modelname is a file of format 

			model_1 	50000	5

		with the first column being the name of the model for a given T and g, second column being
		T and third being log g 


	To find a model with exact values:


	'''

	def __init__(self, lsfilename):
		self.lsfilename = lsfilename
		self.fname = None
		self.t = None
		self.g = None
		self.waves = None
		self.fluxes = None

	def read(self):

		'''
		populates the class instance with the values of T g, and the corresponding 
		names of the model files. Doesn't read in wavelengths and fluxes
		'''

		#self.fname, self.T, self.g

		f, t, g = np.loadtxt(self.lsfilename, unpack=True, dtype = [('name', "|S60"),('t','f8'),('g','f8')])

		self.fname = f
		self.t = t
		self.g = g

		self.waves = []
		self.fluxes = []

		#for i in range(len(f)):
		#	w, fl = np.loadtxt(f[i], unpack=True)

		#	self.waves.append(w)
		#	self.fluxes.append(fl)

		#self.waves = np.array(self.waves)
		#self.fluxes = np.array(self.fluxes)

		return 0


	def find_mod(self, t,g):

		'''
		find a model in the class, return numpy arrays of wavelengths (A)
		and EDDINGTON fluxes (i.e. second moment of radiation field!, an intensity really!)
		'''

		t_find = (self.t == t)
		g_find = (self.g == g)
		find = t_find * g_find

		index = np.where(find == 1)

		ff =  self.fname[index][0]

		w, fl = np.loadtxt(ff, unpack=True)

		return w, 4.0*fl 		# factor of 4 converts from Eddington flux to intensity

	def nn_mod(self, tg_array, mode= "spec"):
		'''
		uses a KDTREE algorithm to find nearest neighbours to array of ts and gs in pairs, shape (2, n)
		'''

		# first construct the tree
		tree = spatial.KDTree(zip(self.t.ravel(), self.g.ravel()))

		# get the locations of nearest neighbours
		results = tree.query(tg_array, k=1)

		if mode == "par":
			# return array with same shape as tg_array containing model values
			return tree.data[results[1]]

		# otherwise, you want to get the actual spectrum for these models
		elif mode == "spec" or mode == "both":
			wavs = []
			fs = []

			for i in range(len(tree.data[results[1]].T[0])):

				t = tree.data[results[1]].T[0][i]
				g = tree.data[results[1]].T[1][i]

				w,f = self.find_mod(t,g) 

				wavs.append( w )
				fs.append(f)

			if mode == "spec":
				return np.array(wavs), np.array(fs)
			elif mode == "both":
				return np.array(wavs), np.array(fs), tree.data[results[1]]

		else:
			print "mode not understood, returning 0"
			return 0




	def near_or_ex(self, tg_array):
		'''
		uses a KDTREE algorithm to find nearest neighbours to array of ts and gs in pairs, shape (2, n).
		Also extrapolate if above highest temperature
		'''


		# first construct the tree
		w, f, dat = self.nn_mod(tg_array, mode= "both")

		tt = (tg_array.T[0])		# temperature is 0th element of array transpose

		tmax = np.max(self.t)

		for i in range(len(tg_array)):
			t = tt[i]

			correction = 1.0

			if t > tmax:

				x1 = np.exp(H_OVER_K * C / (w[i] * ANGSTROM * t)) - 1.0
				x2 = np.exp(H_OVER_K * C / (w[i] * ANGSTROM * tmax)) - 1.0

				correction =  (x2 / x1)


				print "t %8.4e, g %i, correction %8.4e" % (t, tg_array.T[1][i], np.mean(correction))

			f[i] *= correction

		# return array with same shape as tg_array containing model values
		return w,f








'''

def my_interp (par_tree, y, p):
'''
	#	interpolate on an n dimensional grid 
'''

	nparams = len(p)

	nvertices = 2 ** nparams

	weights, locations = tree.query(tg_array, k=nvertices)

	weights = 

	y[results[1]] * results[0] 
	

'''


