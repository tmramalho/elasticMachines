'''
Created on Feb 6, 2013

@author: tiago
'''

import matplotlib.delaunay.triangulate as triang
import numpy as np
from solver import Solver
from automaton import Automaton
import time
import os
from cellSystem import CellSystem
import math
from csolver import dist
import cPickle as cp

class Tissue(object):
	'''
	Main simulation loop manager
	'''

	def __init__(self, k, l, d, nx, ny, ns, ruleID):
		'''
		Constructor
		'''
		self.nx = nx
		self.ny = ny
		self.pt = 0.5*k*math.pow(0.5*d,2)
		self.solver = Solver(k, l, d)
		self.md = 2*d
		self.hd = d/2
		self.ca = Automaton(self.pt)
		self.ca.setRuleID(ruleID)
		self.dt = 0.01
		self.numNewCells = 0
		self.numDead = 0
		self.avEquilib = 0
		self.avSize = 0
		self.numSteps = ns
		self.nIters = 0
		self.dataFolder = os.getcwd() + '/data/em'+str(os.getpid())+str(int(time.time()))+str(ruleID) + "/"
		os.mkdir(self.dataFolder)
		
	def setupDevelopment(self, maxX, maxY):
		self.mx = maxX
		self.my = maxY*4
		cellsX = []
		cellsY = []
		fixed = []
		self.cs = CellSystem()
		for x in xrange(self.nx):
			for y in xrange(self.ny):
				xPos = x/(self.nx-1.0)*maxX + np.random.rand(1)*0.00001
				yPos = y/(self.ny-1.0)*maxY + np.random.rand(1)*0.00001
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
		
	def run(self):
		t1 = time.time()
		for n in range(self.numSteps):
			'''
			Step 1: reach elastic equilibrium
			'''
			if self.numNewCells > 0 or self.numDead > 0:
				self.equilibrate()
			
			'''
			Step 2: save current system state
			'''
			self.saveCellState(n)
			
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
		self.rt = str(t2 - t1)
		self.saveCellState(self.numSteps)
		self.saveParameters()
		
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
		self.avEquilib += i
		self.avSize += self.cs.size
		self.nIters += 1
		
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
			
	def saveParameters(self):
		with open(self.dataFolder+'params.txt', 'w') as f:
			f.write("Run finished at: " + time.strftime("%d %b %Y %H:%M:%S") + "\n")
			f.write("Automaton rule: " + str(self.ca.getRuleID()) + "\n")
			f.write("init num cells x: " + str(self.nx) + "\n")
			f.write("init num cells y: " + str(self.ny) + "\n")
			f.write("pressure threshold: " + str(self.pt) + "\n")
			f.write("lambda: "  + str(self.solver.l) + "\n")
			f.write("stiffness: "  + str(self.solver.k) + "\n")
			f.write("rest length: "  + str(self.solver.d0) + "\n")
			f.write("dt: " + str(self.dt) + "\n")
			f.write("number of ca steps: " + str(self.numSteps) + "\n")
			f.write("max x: " + str(self.mx) + "\n")
			f.write("max y: " + str(self.my) + "\n")
			f.write("average equilibration time: " + str(self.avEquilib/float(self.nIters)) + "\n")
			f.write("average system size: " + str(self.avSize/float(self.nIters)) + "\n")
			f.write("final system size: " + str(self.cs.size) + "\n")
			f.write("runtime in seconds: " + self.rt + "\n")
			
	def saveCellState(self, n):
		with open(self.dataFolder+"step"+str(n).zfill(4)+'.cs', 'wb') as f:
			cp.dump(self.cs, f, protocol=-1)
	