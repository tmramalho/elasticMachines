'''
Created on Feb 6, 2013

@author: tiago
'''

import numpy as np

class Cell(object):
	'''
	Cell properties container
	'''


	def __init__(self, x = 0, y = 0, vx = 0, vy = 0):
		'''
		Constructor
		arrays are such that 0: t+1, 1: t, ... 4: t-3
		this is used in the multistep ode solver
		'''
		self.x = np.array([x for _ in range(5)], dtype=float)
		self.y = np.array([y for _ in range(5)], dtype=float)
		self.vx = np.array([vx for _ in range(5)], dtype=float)
		self.vy = np.array([vy for _ in range(5)], dtype=float)
		self.fx = np.array([0 for _ in range(5)], dtype=float)
		self.fy = np.array([0 for _ in range(5)], dtype=float)
		self.nn = set()
		self.state = 0
		self.newState = 0
		self.fixed = False
		self.V = 0
		
	def addNeighbor(self, n):
		self.nn.add(n)
		
	def forgetNeighbors(self):
		self.nn.clear()
		
	def updatePositions(self):
		self.x = np.roll(self.x, -1)
		self.y = np.roll(self.y, -1)
		self.vx = np.roll(self.vx, -1)
		self.vy = np.roll(self.vy, -1)
		self.fx = np.roll(self.fx, -1)
		self.fy = np.roll(self.fy, -1)
	
	def updateState(self):
		self.state = self.newState