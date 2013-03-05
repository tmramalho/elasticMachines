'''
Created on Feb 7, 2013

@author: tiago
'''

import math
import time
import numpy as np
from csolver import run

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
		
	def runInit(self, cs, dt):
		'''
		Create the initial values for the solver
		using euler method. The first value is left
		as the initial condition
		'''
		for j in xrange(cs.size):
			cs.fx[j, 0], cs.fy[j, 0] = self.calcAcc(cs, j, 0)
		for i in xrange(1, 4):
			for j in xrange(cs.size):
				if cs.fixed[j]:
					continue
				cs.vx[j, i] = cs.vx[j, i-1] + dt*cs.fx[j, i-1]
				cs.vy[j, i] = cs.vy[j, i-1] + dt*cs.fx[j, i-1]
				cs.x[j, i] = cs.x[j, i-1] + dt*cs.vx[j, i-1]
				cs.y[j, i] = cs.y[j, i-1] + dt*cs.vy[j, i-1]
				cs.fx[j, i], cs.fy[j, i] = self.calcAcc(cs, j, i)
	
	def initNewCell(self, cs, dt):
		'''
		Create the previous 3 values for a new cell
		by going back in time
		'''
		j = cs.size - 1
		cs.fx[j, 3], cs.fy[j, 3] = self.calcAcc(cs, j, 3)
		for i in [2, 1, 0]:
			cs.vx[j, i] = cs.vx[j, i+1] - dt*cs.fx[j, i+1]
			cs.vy[j, i] = cs.vy[j, i+1] - dt*cs.fx[j, i+1]
			cs.x[j, i] = cs.x[j, i+1] - dt*cs.vx[j, i+1]
			cs.y[j, i] = cs.y[j, i+1] - dt*cs.vy[j, i+1]
			cs.fx[j, i], cs.fy[j, i] = self.calcAcc(cs, j, i)
		
	def run(self, cs, dt, total):
		'''
		Adams predictor corrector
		'''
		steps = int(total/dt)
		for _ in xrange(steps):
			for j in xrange(cs.size):
				if cs.fixed[j]:
					continue
				cs.vx[j, 4] = cs.vx[j, 3] + dt*(55*cs.fx[j, 3] - 59*cs.fx[j, 2] + 37*cs.fx[j, 1] - 9*cs.fx[j, 0])/24
				cs.vy[j, 4] = cs.vy[j, 3] + dt*(55*cs.fy[j, 3] - 59*cs.fy[j, 2] + 37*cs.fy[j, 1] - 9*cs.fy[j, 0])/24
				cs.x[j, 4]  = cs.x[j, 3]  + dt*(55*cs.vx[j, 3] - 59*cs.vx[j, 2] + 37*cs.vx[j, 1] - 9*cs.vx[j, 0])/24
				cs.y[j, 4]  = cs.y[j, 3]  + dt*(55*cs.vy[j, 3] - 59*cs.vy[j, 2] + 37*cs.vy[j, 1] - 9*cs.vy[j, 0])/24
				cs.fx[j, 4], cs.fy[j, 4] = self.calcAcc(cs, j)
				cs.vx[j, 4] = cs.vx[j, 3] + dt*(9*cs.fx[j, 4] + 19*cs.fx[j, 3] - 5*cs.fx[j, 2] + cs.fx[j, 1])/24
				cs.vy[j, 4] = cs.vy[j, 3] + dt*(9*cs.fy[j, 4] + 19*cs.fy[j, 3] - 5*cs.fy[j, 2] + cs.fy[j, 1])/24
				cs.x[j, 4]  = cs.x[j, 3]  + dt*(9*cs.vx[j, 4] + 19*cs.vx[j, 3] - 5*cs.vx[j, 2] + cs.vx[j, 1])/24
				cs.y[j, 4]  = cs.y[j, 3]  + dt*(9*cs.vy[j, 4] + 19*cs.vy[j, 3] - 5*cs.vy[j, 2] + cs.vy[j, 1])/24
				cs.fx[j, 4], cs.fy[j, 4] = self.calcAcc(cs, j)
			cs.updatePositions()
		
	def calcAcc(self, cs, j, n = 4):
		'''
		calculate acceleration of cell c
		at time t+1 (n=4), t (n=3), 
		        t-1 (n=2), ..., t-3 (n=0)
		'''
		xi = cs.x[j, n]
		yi = cs.y[j, n]
		accx = 0
		accy = 0
		vacc = 0
		for pn in cs.nn[j]:
			xj = cs.x[pn, n]
			yj = cs.y[pn, n]
			dx = xi - xj
			dy = yi - yj
			d = math.sqrt(dx*dx + dy*dy)
			accx += - self.k*dx*(1-self.d0/d)
			accy += - self.k*dy*(1-self.d0/d)
			vacc += 0.5*self.k*math.pow(d - self.d0, 2)
		accx -= self.l*cs.vx[j, n]
		accy -= self.l*cs.vy[j, n]
		cs.V[j] = vacc
		return accx, accy
	
	def fastRun(self, cs, dt, tau):
		nnarr = []
		lenNN = []
		#pack the cell objects onto ndarrays
		for j in xrange(cs.size):
			nnfull = np.zeros(10, dtype = np.int)
			nnset = list(cs.nn[j])
			nnfull[:len(nnset)] = nnset
			nnarr.append(nnfull)
			lenNN.append(len(nnset))
		nnarr = np.array(nnarr, dtype = np.int)
		lenNN = np.array(lenNN, dtype = np.int)
		run(cs.x, cs.y, cs.vx, cs.vy, cs.fx, cs.fy, nnarr, lenNN, cs.V, cs.fixed, dt, tau, self.k, self.l, self.d0)
