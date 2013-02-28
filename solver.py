'''
Created on Feb 7, 2013

@author: tiago
'''

import math
import time

class Solver(object):
	'''
	ODE solver
	'''


	def __init__(self, k, l, d):
		'''
		Constructor
		'''
		self.k = k
		self.l = l
		self.d0 = d
		
	def runInit(self, cells, tCells, dt):
		'''
		Create the initial values for the solver
		using euler method. The first value is left
		as the initial condition
		'''
		for c in tCells:
			c.fx[0], c.fy[0] = self.calcAcc(c, cells, 0)
		for i in xrange(1, 4):
			for c in tCells:
				c.vx[i] = c.vx[i-1] + dt*c.fx[i-1]
				c.vy[i] = c.vy[i-1] + dt*c.fx[i-1]
				c.x[i] = c.x[i-1] + dt*c.vx[i-1]
				c.y[i] = c.y[i-1] + dt*c.vy[i-1]
				c.fx[i], c.fy[i] = self.calcAcc(c, cells, i)
	
	def initNewCell(self, cells, newCell, dt):
		'''
		Create the previous 3 values for a new cell
		by going back in time
		cells only contains old cells so this is safe
		'''
		newCell.fx[3], newCell.fy[3] = self.calcAcc(newCell, cells, 3)
		for i in [2, 1, 0]:
			newCell.vx[i] = newCell.vx[i+1] - dt*newCell.fx[i+1]
			newCell.vy[i] = newCell.vy[i+1] - dt*newCell.fx[i+1]
			newCell.x[i] = newCell.x[i+1] - dt*newCell.vx[i+1]
			newCell.y[i] = newCell.y[i+1] - dt*newCell.vy[i+1]
			newCell.fx[i], newCell.fy[i] = self.calcAcc(newCell, cells, i)
		
	def run(self, cells, tCells, dt, total):
		'''
		Adams predictor corrector
		'''
		t1 = time.time()
		steps = int(total/dt)
		for _ in xrange(steps):
			for c in tCells:
				c.vx[4] = c.vx[3] + dt*(55*c.fx[3] - 59*c.fx[2] + 37*c.fx[1] - 9*c.fx[0])/24
				c.vy[4] = c.vy[3] + dt*(55*c.fy[3] - 59*c.fy[2] + 37*c.fy[1] - 9*c.fy[0])/24
				c.x[4]  = c.x[3]  + dt*(55*c.vx[3] - 59*c.vx[2] + 37*c.vx[1] - 9*c.vx[0])/24
				c.y[4]  = c.y[3]  + dt*(55*c.vy[3] - 59*c.vy[2] + 37*c.vy[1] - 9*c.vy[0])/24
				c.fx[4], c.fy[4] = self.calcAcc(c, cells)
				c.vx[4] = c.vx[3] + dt*(9*c.fx[4] + 19*c.fx[3] - 5*c.fx[2] + c.fx[1])/24
				c.vy[4] = c.vy[3] + dt*(9*c.fy[4] + 19*c.fy[3] - 5*c.fy[2] + c.fy[1])/24
				c.x[4]  = c.x[3]  + dt*(9*c.vx[4] + 19*c.vx[3] - 5*c.vx[2] + c.vx[1])/24
				c.y[4]  = c.y[3]  + dt*(9*c.vy[4] + 19*c.vy[3] - 5*c.vy[2] + c.vy[1])/24
				c.fx[4], c.fy[4] = self.calcAcc(c, cells)
			for c in tCells:
				c.updatePositions()
		t2 = time.time()
		print("Time = " + str(t2 - t1))
		
	def calcAcc(self, c, cells, n = 4):
		'''
		calculate acceleration of cell c
		at time t+1 (n=4), t (n=3), 
		        t-1 (n=2), ..., t-3 (n=0)
		'''
		xi = c.x[n]
		yi = c.y[n]
		accx = 0
		accy = 0
		c.V = 0
		for pn in c.nn:
			xj = cells[pn].x[n]
			yj = cells[pn].y[n]
			dx = xi - xj
			dy = yi - yj
			d = math.sqrt(dx*dx + dy*dy)
			accx += - self.k*dx*(1-self.d0/d)
			accy += - self.k*dy*(1-self.d0/d)
			c.V += 0.5*self.k*math.pow(d - self.d0, 2)
		accx -= self.l*c.vx[n]
		accy -= self.l*c.vy[n]
		return accx, accy