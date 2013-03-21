'''
Created on Mar 3, 2013

@author: tiago
'''

import numpy as np

class CellSystem():
	'''
	A container for all the cells
	'''


	def __init__(self):
		'''
		Setup properties
		'''
		
	def initializeCells(self, xPos, yPos, fixed):
		'''
		arrays are such that 0: t+1, 1: t, ... 4: t-3
		this is used in the multistep ode solver
		'''
		numCells = len(xPos)
		self.x = np.reshape(np.repeat(xPos, 5), (numCells, 5))
		self.y = np.reshape(np.repeat(yPos, 5), (numCells, 5))
		self.vx = np.zeros((numCells, 5), dtype=float)
		self.vy = np.zeros((numCells, 5), dtype=float)
		self.fx = np.zeros((numCells, 5), dtype=float)
		self.fy = np.zeros((numCells, 5), dtype=float)
		self.nn = [set() for _ in xrange(numCells)]
		self.state = np.zeros(numCells, dtype=int)
		self.newState = np.zeros(numCells, dtype=int)
		self.fixed = np.array(fixed, dtype=int)
		self.V = np.zeros(numCells, dtype=float)
		self.size = self.x.shape[0]
		
	def addCell(self, xPos, yPos, fixed):
		self.x = np.concatenate([self.x, [np.repeat(xPos, 5)]])
		self.y = np.concatenate([self.y, [np.repeat(yPos, 5)]])
		self.vx = np.concatenate([self.vx, [np.zeros(5)]])
		self.vy = np.concatenate([self.vy, [np.zeros(5)]])
		self.fx = np.concatenate([self.fx, [np.zeros(5)]])
		self.fy = np.concatenate([self.fy, [np.zeros(5)]])
		self.nn.append(set())
		self.state = np.concatenate([self.state, [0]])
		self.newState = np.concatenate([self.newState, [0]])
		self.fixed = np.concatenate([self.fixed, [False]])
		self.V = np.concatenate([self.V, [0]])
		self.size = self.x.shape[0]
		
	def deleteCells(self, targets):
		self.x = np.delete(self.x, targets, axis=0)
		self.y = np.delete(self.y, targets, axis=0)
		self.vx = np.delete(self.vx, targets, axis=0)
		self.vy = np.delete(self.vy, targets, axis=0)
		self.fx = np.delete(self.fx, targets, axis=0)
		self.fy = np.delete(self.fy, targets, axis=0)
		self.state = np.delete(self.state, targets, axis=0)
		self.newState = np.delete(self.newState, targets, axis=0)
		self.fixed = np.delete(self.fixed, targets, axis=0)
		self.V = np.delete(self.V, targets, axis=0)
		for t in sorted(targets, reverse=True):
			self.nn.pop(t)
		self.size = self.x.shape[0]
		
	def deleteRogueCells(self, maxY):
		tol = 0.001
		bx = 0 - tol
		tx = 1 + tol
		by = 0 - tol
		ty = maxY + tol
		numCells = self.x.shape[0]
		goners = set()
		for j in xrange(numCells):
			if self.x[j, 3] < bx or self.x[j, 3] > tx or self.y[j, 3] < by or self.y[j, 3] > ty:
				goners.add(j)
		self.deleteCells(list(goners))
		self.size = self.x.shape[0]
		return len(goners)
		
	def updatePositions(self):
		self.x = np.roll(self.x, -1, axis=1)
		self.y = np.roll(self.y, -1, axis=1)
		self.vx = np.roll(self.vx, -1, axis=1)
		self.vy = np.roll(self.vy, -1, axis=1)
		self.fx = np.roll(self.fx, -1, axis=1)
		self.fy = np.roll(self.fy, -1, axis=1)
	
	def addNeighbor(self, j, n):
		self.nn[j].add(n)
	
	def forgetAllNeighbors(self):
		for nSet in self.nn:
			nSet.clear()
		
	def forgetNeighbors(self, j):
		self.nn[j].clear()
	
	def updateState(self):
		#swap the two arrays
		self.state, self.newState = self.newState, self.state
		