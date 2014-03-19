from disky_const import *
from scipy import integrate, spatial
import numpy as np


class model:

	'''
	A class to store a set of models
	'''

	def __init__(self, lsfilename):
		self.lsfilename = lsfilename
		self.fname = None
		self.t = None
		self.g = None
		self.waves = None
		self.fluxes = None

	def read(self):

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

		t_find = (self.t == t)
		g_find = (self.g == g)
		find = t_find * g_find

		index = np.where(find == 1)

		ff =  self.fname[index][0]

		w, fl = np.loadtxt(ff, unpack=True)

		return w, fl

	def nearrest_mod(self, tg_array):
		'''
		uses a KDTREE algorithm to find nearest neighbours to array of ts and gs in pairs, shape (2, n)
		'''

		# first construct the tree
		tree = spatial.KDTree(zip(t.ravel(), g.ravel()))

		# get the locations of nearest neighbours
		results = tree.query(tg_array, k=1)

		# return array with same shape as tg_array containing model values
		return tree.data[results[1]]









