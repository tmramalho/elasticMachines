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
		self.transition = np.zeros(8)
		self.transition[:] = [0,1,1,0,1,1,1,0]
		self.growth = np.zeros(8)
		self.growth[:] = [0,0,1,1,0,0,1,0]
		self.pt = pt
		
	def genpos(self, s, x, p):
		'''
		generate position in transition table
		corresponding to inputs s,x,p
		s int [0,1]
		x int [0,1]
		p int [0,1]
		'''
		return s<<2 | x<<1 | p
	
	def getRuleID(self):
		rid = 0
		for i in xrange(8):
			rid += self.transition[i]*np.power(2,i)
			rid += self.growth[i]*np.power(2,i+8)
		return int(rid)
	
	def setRuleID(self, rid):
		bits = np.binary_repr(rid, width=16)[::-1]
		print bits
		for i in xrange(8):
			self.transition[i] = bits[i]
			self.growth[i] = bits[i+8]
	
	def evolve(self, cs):
		self.newCells = set()
		for j in xrange(cs.size):
			acc = 0
			for pn in cs.nn[j]:
				acc += cs.state[pn]
			s = 0 if acc < 3 else 1
			p = 0 if cs.V[j] < self.pt else 1
			x = cs.state[j]
			i = self.genpos(s, x, p)
			try:
				cs.newState[j] = self.transition[i]
			except IndexError:
				print "whoopsie", i, p, x, acc, len(cs.nn[j])
				exit()
			if self.growth[i] == 1:
				self.newCells.add((cs.x[j, 3], cs.y[j, 3]))
		cs.updateState()
		