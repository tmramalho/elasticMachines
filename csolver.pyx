from __future__ import division
from libc.math cimport sqrt
from libc.math cimport pow
from libc.stdio cimport printf
import numpy as np
cimport numpy as np
cimport cython

FLDT = np.float
ctypedef np.float_t FLDT_t
INTDT = np.int
ctypedef np.int_t INTDT_t

@cython.boundscheck(False)
@cython.cdivision(True)
def run(
	np.ndarray[FLDT_t, ndim=2] xarr, 
	np.ndarray[FLDT_t, ndim=2] yarr, 
	np.ndarray[FLDT_t, ndim=2] vxarr, 
	np.ndarray[FLDT_t, ndim=2] vyarr, 
	np.ndarray[FLDT_t, ndim=2] fxarr, 
	np.ndarray[FLDT_t, ndim=2] fyarr, 
	np.ndarray[INTDT_t, ndim=2] nnarr, 
	np.ndarray[INTDT_t, ndim=1] lenNN,
	np.ndarray[FLDT_t, ndim=1] varr, 
	np.ndarray[INTDT_t, ndim=1] fiarr, 
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
	cdef unsigned int four = 4
	cdef unsigned int three = 3
	cdef unsigned int two = 2
	cdef unsigned int one = 1
	cdef unsigned int zero = 0
	for i in range(steps):
		if fiarr[j] == 1:
			continue
		else:
			printf("wazoo")
		for j in range(nCells):
			vxarr[j][four] = vxarr[j][three] + dt*(55*fxarr[j][three] - 59*fxarr[j][two] + 37*fxarr[j][one] - 9*fxarr[j][zero])/24
			vyarr[j][four] = vyarr[j][three] + dt*(55*fyarr[j][three] - 59*fyarr[j][two] + 37*fyarr[j][one] - 9*fyarr[j][zero])/24
			xarr[j][four]  = xarr[j][three]  + dt*(55*vxarr[j][three] - 59*vxarr[j][two] + 37*vxarr[j][one] - 9*vxarr[j][zero])/24
			yarr[j][four]  = yarr[j][three]  + dt*(55*vyarr[j][three] - 59*vyarr[j][two] + 37*vyarr[j][one] - 9*vyarr[j][zero])/24
			#begin inline calcAcc
			acc[zero] = 0
			acc[one] = 0
			xi = xarr[j][four]
			yi = yarr[j][four]
			vacc = 0
			for n in range(lenNN[j]):
				pn = nnarr[j][n]
				xn = xarr[pn][four]
				yn = yarr[pn][four]
				dx = xi - xn
				dy = yi - yn
				d = sqrt(dx*dx + dy*dy)
				acc[zero] -= k*dx*(1-d0/d)
				acc[one] -= k*dy*(1-d0/d)
				vacc += 0.5*k*pow(d - d0, 2)
			varr[j] = vacc
			acc[zero] -= l*vxarr[j][four]
			acc[one] -= l*vyarr[j][four]
			fxarr[j][four] = acc[zero]
			fyarr[j][four] = acc[one]
			#end inline calcAcc
			vxarr[j][four] = vxarr[j][three] + dt*(9*fxarr[j][four] + 19*fxarr[j][three] - 5*fxarr[j][two] + fxarr[j][one])/24
			vyarr[j][four] = vyarr[j][three] + dt*(9*fyarr[j][four] + 19*fyarr[j][three] - 5*fyarr[j][two] + fyarr[j][one])/24
			xarr[j][four]  = xarr[j][three]  + dt*(9*vxarr[j][four] + 19*vxarr[j][three] - 5*vxarr[j][two] + vxarr[j][one])/24
			yarr[j][four]  = yarr[j][three]  + dt*(9*vyarr[j][four] + 19*vyarr[j][three] - 5*vyarr[j][two] + vyarr[j][one])/24
			#begin inline calcAcc
			acc[zero] = 0
			acc[one] = 0
			xi = xarr[j][four]
			yi = yarr[j][four]
			vacc = 0
			for n in range(lenNN[j]):
				pn = nnarr[j][n]
				xn = xarr[pn][four]
				yn = yarr[pn][four]
				dx = xi - xn
				dy = yi - yn
				d = sqrt(dx*dx + dy*dy)
				acc[zero] -= k*dx*(1-d0/d)
				acc[one] -= k*dy*(1-d0/d)
				vacc += 0.5*k*pow(d - d0, 2)
			varr[j] = vacc
			acc[zero] -= l*vxarr[j][four]
			acc[one] -= l*vyarr[j][four]
			fxarr[j][four] = acc[zero]
			fyarr[j][four] = acc[one]
			#end inline calcAcc
		for j in range(nCells):
			for n in range(4):
				xarr[j][n] = xarr[j][n+one]
				yarr[j][n] = yarr[j][n+one]
				vxarr[j][n] = vxarr[j][n+one]
				vyarr[j][n] = vyarr[j][n+one]
				fxarr[j][n] = fxarr[j][n+one]
				fyarr[j][n] = fyarr[j][n+one]

@cython.boundscheck(False)
@cython.cdivision(True)
cdef inline double *calcAcc(
	double** xarr, 
	double** yarr, 
	double** vxarr, 
	double** vyarr, 
	int** nnarr,
	int* lenNN, 
	double* varr, 
	unsigned int j, double k, double l, double d0):
	'''
	calculate acceleration of cell c
	at time t+1 (n=4), t (n=3), 
	        t-1 (n=2), ..., t-3 (n=0)
	'''
	cdef unsigned int zero = 0
	cdef unsigned int one = 1
	cdef unsigned int four = 4
	cdef double acc[2]
	cdef double xi = xarr[j][four]
	cdef double yi = yarr[j][four]
	cdef double vacc = 0
	cdef double xn
	cdef double yn
	cdef double dx
	cdef double dy
	cdef double d
	cdef unsigned int n
	cdef unsigned int pn
	acc[zero] = 0
	acc[one] = 0
	for n in range(lenNN[j]):
		pn = nnarr[j][n]
		xn = xarr[pn][four]
		yn = yarr[pn][four]
		dx = xi - xn
		dy = yi - yn
		d = sqrt(dx*dx + dy*dy)
		acc[zero] -= k*dx*(1-d0/d)
		acc[one] -= k*dy*(1-d0/d)
		vacc += 0.5*k*pow(d - d0, 2)
	varr[j] = vacc
	acc[zero] -= l*vxarr[j][four]
	acc[one] -= l*vyarr[j][four]
	return acc