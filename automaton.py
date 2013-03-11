'''
Created on Feb 11, 2013

@author: tiago
'''

import numpy as np

class Automaton(object):
	'''
	classdocs
	'''


	def __init__(self, pt):
		'''
		Constructor initializes default values
		for automaton simulation
		'''
		self.pt = 1 #pressure threshold
		self.transition = np.zeros(32)
		self.transition[self.genpos(1,0,0)] = 1
		self.transition[self.genpos(1,0,1)] = 1
		self.transition[self.genpos(4,0,0)] = 1
		self.transition[self.genpos(4,0,1)] = 1
		self.transition[self.genpos(6,0,0):] = 1
		self.growth = np.zeros(32)
		self.growth[self.genpos(4,0,0)] = 1
		self.growth[self.genpos(4,0,1)] = 1
		self.pt = pt
		
	def genpos(self, s, x, p):
		'''
		generate position in transition table
		corresponding to inputs s,x,p
		s int [0,7]
		x int [0,1]
		p int [0,1]
		'''
		return s<<2 | x<<1 | p
	
	def evolve(self, cs):
		self.newCells = set()
		for j in xrange(cs.size):
			acc = 0
			for pn in cs.nn[j]:
				acc += cs.state[pn]
			acc = acc if acc < 8 else 7
			p = 0 if cs.V[j] < self.pt else 1
			x = cs.state[j]
			i = self.genpos(acc, x, p)
			try:
				cs.newState[j] = self.transition[i]
			except IndexError:
				print "whoopsie", i, p, x, acc, len(cs.nn[j])
				exit()
			if self.growth[i] == 1:
				self.newCells.add((cs.x[j, 3], cs.y[j, 3]))
		cs.updateState()
		