from __future__ import division
from libc.math cimport sqrt
from libc.math cimport pow
from libc.stdio cimport printf
import numpy as np
cimport numpy as np
cimport cython

FLDT = np.float
INTDT = np.int
ctypedef np.float_t FLDT_t
ctypedef np.int_t INTDT_t

@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def run(
	np.ndarray[FLDT_t, ndim=2] xarr not None, 
	np.ndarray[FLDT_t, ndim=2] yarr not None, 
	np.ndarray[FLDT_t, ndim=2] vxarr not None, 
	np.ndarray[FLDT_t, ndim=2] vyarr not None, 
	np.ndarray[FLDT_t, ndim=2] fxarr not None, 
	np.ndarray[FLDT_t, ndim=2] fyarr not None, 
	np.ndarray[INTDT_t, ndim=2] nnarr not None, 
	np.ndarray[INTDT_t, ndim=1] lenNN not None,
	np.ndarray[FLDT_t, ndim=1] varr not None, 
	np.ndarray[INTDT_t, ndim=1] fiarr not None, 
	double dt, double total, double k, double l, double d0):
	'''
	Adams predictor corrector
	'''
	cdef int steps = int(total/dt)
	cdef double * f
	cdef int nCells = xarr.shape[0]
	cdef unsigned int i,j,n,numnn,pn
	cdef double xi,xj,yi,yj,vacc,dx,dy,d
	cdef double acc[2]
	cdef double test = 0
	for i in range(steps):
		for j in range(nCells):
			if fiarr[j] == 1:
				continue
			test = vxarr[j, 3]/2
			vxarr[j, 4] = vxarr[j, 3] + dt*(55*fxarr[j, 3] - 59*fxarr[j, 2] + 37*fxarr[j, 1] - 9*fxarr[j, 0])/24
			vyarr[j, 4] = vyarr[j, 3] + dt*(55*fyarr[j, 3] - 59*fyarr[j, 2] + 37*fyarr[j, 1] - 9*fyarr[j, 0])/24
			xarr[j, 4]  = xarr[j, 3]  + dt*(55*vxarr[j, 3] - 59*vxarr[j, 2] + 37*vxarr[j, 1] - 9*vxarr[j, 0])/24
			yarr[j, 4]  = yarr[j, 3]  + dt*(55*vyarr[j, 3] - 59*vyarr[j, 2] + 37*vyarr[j, 1] - 9*vyarr[j, 0])/24
			#begin inline calcAcc
			acc[0] = 0
			acc[1] = 0
			xi = xarr[j, 4]
			yi = yarr[j, 4]
			vacc = 0
			for n in range(lenNN[j]):
				pn = nnarr[j, n]
				xn = xarr[pn, 4]
				yn = yarr[pn, 4]
				dx = xi - xn
				dy = yi - yn
				d = sqrt(dx*dx + dy*dy)
				acc[0] -= k*dx*(1-d0/d)
				acc[1] -= k*dy*(1-d0/d)
				vacc += 0.5*k*pow(d - d0, 2)
			varr[j] = vacc
			acc[0] -= l*vxarr[j, 4]
			acc[1] -= l*vyarr[j, 4]
			fxarr[j, 4] = acc[0]
			fyarr[j, 4] = acc[1]
			#end inline calcAcc
			vxarr[j, 4] = vxarr[j, 3] + dt*(9*fxarr[j, 4] + 19*fxarr[j, 3] - 5*fxarr[j, 2] + fxarr[j, 1])/24
			vyarr[j, 4] = vyarr[j, 3] + dt*(9*fyarr[j, 4] + 19*fyarr[j, 3] - 5*fyarr[j, 2] + fyarr[j, 1])/24
			xarr[j, 4]  = xarr[j, 3]  + dt*(9*vxarr[j, 4] + 19*vxarr[j, 3] - 5*vxarr[j, 2] + vxarr[j, 1])/24
			yarr[j, 4]  = yarr[j, 3]  + dt*(9*vyarr[j, 4] + 19*vyarr[j, 3] - 5*vyarr[j, 2] + vyarr[j, 1])/24
			#begin inline calcAcc
			acc[0] = 0
			acc[1] = 0
			xi = xarr[j, 4]
			yi = yarr[j, 4]
			vacc = 0
			for n in range(lenNN[j]):
				pn = nnarr[j, n]
				xn = xarr[pn, 4]
				yn = yarr[pn, 4]
				dx = xi - xn
				dy = yi - yn
				d = sqrt(dx*dx + dy*dy)
				acc[0] -= k*dx*(1-d0/d)
				acc[1] -= k*dy*(1-d0/d)
				vacc += 0.5*k*pow(d - d0, 2)
			varr[j] = vacc
			acc[0] -= l*vxarr[j, 4]
			acc[1] -= l*vyarr[j, 4]
			fxarr[j, 4] = acc[0]
			fyarr[j, 4] = acc[1]
			#end inline calcAcc
		for j in range(nCells):
			for n in range(4):
				xarr[j, n] = xarr[j, n+1]
				yarr[j, n] = yarr[j, n+1]
				vxarr[j, n] = vxarr[j, n+1]
				vyarr[j, n] = vyarr[j, n+1]
				fxarr[j, n] = fxarr[j, n+1]
				fyarr[j, n] = fyarr[j, n+1]
