'''
Created on Feb 11, 2013

@author: tiago
'''

class Automaton(object):
	'''
	classdocs
	'''


	def __init__(self):
		'''
		Constructor
		'''
	
	def stateEvolve(self, cs):
		self.newCells = set()
		for j in xrange(cs.size):
			acc = 0
			for pn in cs.nn[j]:
				acc += cs.state[pn]
			if acc == 1:
				if cs.state[j] == 1:
					cs.newState[j] = 0
				else:
					cs.newState[j] = 1
			elif acc > 3 and acc < 5:
				if cs.state[j] == 0:
					if not cs.fixed[j]:
						self.newCells.add((cs.x[j, 3], cs.y[j, 3]))
					cs.newState[j] = 1
				else:
					cs.newState[j] = 0
			elif acc > 5:
				cs.newState[j] = 0
		cs.updateState()