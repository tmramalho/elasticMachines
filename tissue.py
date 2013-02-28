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
		self.fCells = []
		self.tCells = []
		for x in xrange(self.nx):
			for y in xrange(self.ny):
				xPos = x/(self.nx-1.0)*maxX
				yPos = y/(self.ny-1.0)*maxY
				c = Cell(xPos, yPos)
				#or y == self.ny-1
				if(x == 0 or y == 0 or x == self.nx-1):
					c.fixed = True
					self.fCells.append(c)
				else:
					c.x += (np.random.rand(1)-0.5)*self.hd
					c.y += (np.random.rand(1)-0.5)*self.hd
					self.tCells.append(c)
		self.cells = self.fCells + self.tCells
		self.calcDelaunay()
		self.solver.runInit(self.cells, self.tCells, 0.01)
		self.tCells[-3].state = 1
		self.plotter.plotCellStates(self.cells, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes)
		self.plotter.next()
		
	def run(self):
		for _ in range(20):
			'''
			Step 1a: evolve until equilibrium
			'''
			self.solver.run(self.cells, self.tCells, self.dt, 2)
			
			'''
			Step 1b: rebuild links
			'''
			self.calcDelaunay()
			
			'''
			Step 2: plot stuff
			'''
			#self.plotter.plotDelaunay(self.xp, self.yp, self.tt.triangle_nodes)
			#self.plotter.plotVoronoi(self.cells, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes, self.tt.triangle_nodes)
			self.plotter.plotCellStates(self.cells, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes)
			self.plotter.drawCells(self.cells, self.tt.triangle_nodes, self.xp, self.yp)
			self.plotter.next()
			
			'''
			Step 3: update state and add cells
			'''
			self.ca.stateEvolve(self.cells)
			
			'''
			Step 4: add cells (tCells)
			'''
			self.grow()
			
			'''
			Step 5: delete cells which went outside the tissue (tCells)
			'''
			self.deleteRogueCells()
			
			'''
			Step 6: update array with all cells
			'''
			self.cells = self.fCells + self.tCells
			
			'''
			Step 7: rebuild links
			'''
			self.calcDelaunay()
		
	def calcDelaunay(self):
		#delaunay
		xp = []
		yp = []
		for cell in self.cells:
			cell.forgetNeighbors()
			xp.append(cell.x[3])
			yp.append(cell.y[3])
		self.xp = np.array(xp)
		self.yp = np.array(yp)
		self.tt = triang.Triangulation(xp,yp)
		self.edges = [e for e in self.tt.edge_db if self.dist(e,xp,yp) < self.md]
		for e in self.edges:
			self.cells[e[0]].addNeighbor(e[1])
			self.cells[e[1]].addNeighbor(e[0])
		
	def dist(self, e, xp, yp):
		dx = xp[e[0]] - xp[e[1]]
		dy = yp[e[0]] - yp[e[1]]
		return np.sqrt(dx*dx+dy*dy)
	
	def grow(self):
		for nc in self.ca.newCells:
			r = self.hd*(np.random.rand(2) - 0.5)
			c = Cell(nc[0]+r[0], nc[1]+r[1])
			self.tCells.append(c)
			self.solver.initNewCell(self.cells, c, self.dt)
			
	def deleteRogueCells(self):
		self.tCells[:] = [c for c in self.tCells if c.x[3] >= 0 and c.x[3] <= 1 and c.y[3] >= 0]
