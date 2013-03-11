'''
Created on Feb 6, 2013

@author: tiago
'''

import matplotlib.delaunay.triangulate as triang
from plotter import Plotter
import numpy as np
from solver import Solver
from automaton import Automaton
import time
from cellSystem import CellSystem
import math
from csolver import dist

class Tissue(object):
	'''
	Main simulation loop manager
	'''


	def __init__(self, k, l, d, nx, ny, plot=False):
		'''
		Constructor
		'''
		self.nx = nx
		self.ny = ny
		self.pt = 0.5*k*math.pow(0.5*d,2)
		self.plotter = Plotter(d, self.pt)
		self.solver = Solver(k, l, d)
		self.md = 2*d
		self.hd = d/2
		self.ca = Automaton(self.pt)
		self.dt = 0.01
		self.plotAllSteps = plot
		self.numNewCells = 0
		self.numDead = 0
		
	def setupDevelopment(self, maxX, maxY):
		self.plotter.setLims(maxX, maxY)
		self.my = maxY*4
		print "end:", self.my
		cellsX = []
		cellsY = []
		fixed = []
		self.cs = CellSystem()
		for x in xrange(self.nx):
			for y in xrange(self.ny):
				xPos = x/(self.nx-1.0)*maxX
				yPos = y/(self.ny-1.0)*maxY
				#or y == self.ny-1
				if(x == 0 or y == 0 or x == self.nx-1):
					fixed.append(True)
				else:
					#xPos += (np.random.rand(1)-0.5)*self.hd
					#yPos += (np.random.rand(1)-0.5)*self.hd
					fixed.append(False)
				cellsX.append(xPos)
				cellsY.append(yPos)
		self.cs.initializeCells(cellsX, cellsY, fixed)
		self.numNewCells = self.nx*self.ny
		self.calcDelaunay()
		self.solver.runInit(self.cs, 0.01)
		self.cs.state[140] = 1
		self.plotter.plotCellStates(self.cs, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes)
		self.plotter.next()
		
	def run(self):
		t1 = time.time()
		for _ in range(50):
			if self.numNewCells > 0 or self.numDead > 0:
				self.equilibrate()
			
			if self.plotAllSteps:
				'''
				Step 2: plot stuff
				'''
				self.plotter.plotCellStates(self.cs, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes)
				self.plotter.drawCells(self.cs, self.tt.triangle_nodes, self.xp, self.yp)
				self.plotter.next()
			
			'''
			Step 3: update state and add cells
			'''
			self.ca.evolve(self.cs)
			self.numNewCells = len(self.ca.newCells)
			
			'''
			Step 4: add cells (tCells)
			'''
			self.grow()
			
			'''
			Step 5: delete cells which went outside the tissue (tCells)
			'''
			self.numDead = self.cs.deleteRogueCells(self.my)
		
		self.equilibrate()	
		t2 = time.time()
		print("Time 1a = " + str(t2 - t1))
		self.plotter.plotCellStates(self.cs, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes)
		self.plotter.drawCells(self.cs, self.tt.triangle_nodes, self.xp, self.yp)
		self.plotter.next()
		
	def equilibrate(self):
		i = 0
		while True:
			self.calcDelaunay()
			#self.solver.run(self.cs, self.dt, 2)
			self.solver.fastRun(self.cs, self.dt, 2)
			tx = np.allclose(self.cs.x[:,3], self.xp)
			ty = np.allclose(self.cs.y[:,3], self.yp)
			if ty and tx:
				break
			i += 1
		print "equilib:", i, self.cs.size
		
	def calcDelaunay(self):
		#delaunay
		self.xp = np.copy(self.cs.x[:,3])
		self.yp = np.copy(self.cs.y[:,3])
		self.cs.forgetAllNeighbors()
		self.tt = triang.Triangulation(self.xp,self.yp)
		self.edges = [e for e in self.tt.edge_db if dist(e,self.xp,self.yp) < self.md]
		self.numEdges = len(self.edges)
		for e in self.edges:
			self.cs.addNeighbor(e[0], e[1])
			self.cs.addNeighbor(e[1], e[0])
	
	def grow(self):
		for nc in self.ca.newCells:
			r = self.hd*(np.random.rand(2) - 0.5)
			self.cs.addCell(nc[0]+r[0], nc[1]+r[1], False)
			self.solver.initNewCell(self.cs, self.dt)
	