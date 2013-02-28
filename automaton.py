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
	
	def stateEvolve(self, cells):
		self.newCells = set()
		for c in cells:
			acc = 0
			for pn in c.nn:
				acc += cells[pn].state
			if acc == 1:
				if c.state == 1:
					c.newState = 0
				else:
					c.newState = 1
			elif acc > 3 and acc < 5:
				if c.state == 0:
					if not c.fixed:
						self.newCells.add((c.x[3], c.y[3]))
					c.newState = 1
				else:
					c.newState = 0
			elif acc > 5:
				c.newState = 0
		for c in cells:
			c.updateState()