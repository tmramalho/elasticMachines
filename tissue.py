'''
Created on Feb 6, 2013

@author: tiago
'''

from cell import Cell
import matplotlib.delaunay.triangulate as triang
from plotter import Plotter
import numpy as np
from solver import Solver
from automaton import Automaton
import time
from cellSystem import CellSystem

class Tissue(object):
	'''
	Main simulation loop manager
	'''


	def __init__(self, k, l, d, nx, ny):
		'''
		Constructor
		'''
		self.nx = nx
		self.ny = ny
		self.plotter = Plotter(d)
		self.solver = Solver(k, l, d)
		self.md = 2*d
		self.hd = d/2
		self.ca = Automaton()
		self.dt = 0.01
		
	def setupDevelopment(self, maxX, maxY):
		self.plotter.setLims(maxX, maxY)
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
		self.calcDelaunay()
		self.solver.runInit(self.cs, 0.01)
		self.cs.state[22] = 1
		self.plotter.plotCellStates(self.cs, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes)
		self.plotter.next()
		
	def run(self):
		for _ in range(20):
			'''
			Step 1a: evolve until equilibrium
			'''
			t1 = time.time()
			#self.solver.run(self.cs, self.dt, 2)
			self.solver.fastRun(self.cs, self.dt, 2)
			t2 = time.time()
			print("Time 1a = " + str(t2 - t1))
			
			'''
			Step 1b: rebuild links
			'''
			self.calcDelaunay()
			
			'''
			Step 2: plot stuff
			'''
			#self.plotter.plotDelaunay(self.xp, self.yp, self.tt.triangle_nodes)
			#self.plotter.plotVoronoi(self.cs, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes, self.tt.triangle_nodes)
			self.plotter.plotCellStates(self.cs, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes)
			self.plotter.drawCells(self.cs, self.tt.triangle_nodes, self.xp, self.yp)
			self.plotter.next()
			
			'''
			Step 3: update state and add cells
			'''
			self.ca.stateEvolve(self.cs)
			
			'''
			Step 4: add cells (tCells)
			'''
			self.grow()
			
			'''
			Step 5: delete cells which went outside the tissue (tCells)
			'''
			self.cs.deleteRogueCells()
			
			'''
			Step 7: rebuild links
			'''
			self.calcDelaunay()
		
	def calcDelaunay(self):
		#delaunay
		self.xp = self.cs.x[:,3]
		self.yp = self.cs.y[:,3]
		self.cs.forgetAllNeighbors()
		self.tt = triang.Triangulation(self.xp,self.yp)
		self.edges = [e for e in self.tt.edge_db if self.dist(e,self.xp,self.yp) < self.md]
		for e in self.edges:
			self.cs.addNeighbor(e[0], e[1])
			self.cs.addNeighbor(e[1], e[0])
		
	def dist(self, e, xp, yp):
		dx = xp[e[0]] - xp[e[1]]
		dy = yp[e[0]] - yp[e[1]]
		return np.sqrt(dx*dx+dy*dy)
	
	def grow(self):
		for nc in self.ca.newCells:
			r = self.hd*(np.random.rand(2) - 0.5)
			self.cs.addCell(nc[0]+r[0], nc[1]+r[1], False)
			self.solver.initNewCell(self.cs, self.dt)
	