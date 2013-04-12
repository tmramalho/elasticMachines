/*
 * Solver.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: tiago
 */

#include "Solver.h"

void Solver::run(double dt, double total, double k, double l, double d0,
		std::vector<double>& xarr, std::vector<double>& yarr,
		std::vector<double>& vxarr, std::vector<double>& vyarr,
		std::vector<double>& fxarr, std::vector<double>& fyarr,
		std::vector<std::set<int>*>& nnarr, std::vector<double>& varr,
		std::vector<bool>& fiarr) {

	int steps = (int) ((double)total/dt);
	int nCells = fiarr.size();
	double vacc, xi, yi, d, dx, dy, xn, yn;
	int pn;
	for (int i = 0; i < steps; i++) {
		for (int j = 0; j < nCells; j++) {
			if (fiarr[j]) continue;
			vxarr[j*4 +4] = vxarr[j*4 +3] + dt*(55*fxarr[j*4 +3] - 59*fxarr[j*4 +2] + 37*fxarr[j*4 +1] - 9*fxarr[j*4 +0])/24;
			vyarr[j*4 +4] = vyarr[j*4 +3] + dt*(55*fyarr[j*4 +3] - 59*fyarr[j*4 +2] + 37*fyarr[j*4 +1] - 9*fyarr[j*4 +0])/24;
			xarr[j*4 +4]  = xarr[j*4 +3]  + dt*(55*vxarr[j*4 +3] - 59*vxarr[j*4 +2] + 37*vxarr[j*4 +1] - 9*vxarr[j*4 +0])/24;
			yarr[j*4 +4]  = yarr[j*4 +3]  + dt*(55*vyarr[j*4 +3] - 59*vyarr[j*4 +2] + 37*vyarr[j*4 +1] - 9*vyarr[j*4 +0])/24;
			calcAcc(4, 3, j, dt, k, l, d0, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
			vxarr[j*4 +4] = vxarr[j*4 +3] + dt*(9*fxarr[j*4 +4] + 19*fxarr[j*4 +3] - 5*fxarr[j*4 +2] + fxarr[j*4 +1])/24;
			vyarr[j*4 +4] = vyarr[j*4 +3] + dt*(9*fyarr[j*4 +4] + 19*fyarr[j*4 +3] - 5*fyarr[j*4 +2] + fyarr[j*4 +1])/24;
			xarr[j*4 +4]  = xarr[j*4 +3]  + dt*(9*vxarr[j*4 +4] + 19*vxarr[j*4 +3] - 5*vxarr[j*4 +2] + vxarr[j*4 +1])/24;
			yarr[j*4 +4]  = yarr[j*4 +3]  + dt*(9*vyarr[j*4 +4] + 19*vyarr[j*4 +3] - 5*vyarr[j*4 +2] + vyarr[j*4 +1])/24;
			calcAcc(4, 3, j, dt, k, l, d0, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
		}
		for (int j = 0; j < nCells; j++) {
			if (fiarr[j]) continue;
			for (int n = 0; n < 4; n++) {
				//shift everyone back one
				xarr[j*4 +n] = xarr[j*4 + n+1];
				yarr[j*4 +n] = yarr[j*4 + n+1];
				vxarr[j*4 +n] = vxarr[j*4 + n+1];
				vyarr[j*4 +n] = vyarr[j*4 + n+1];
				fxarr[j*4 +n] = fxarr[j*4 + n+1];
				fyarr[j*4 +n] = fyarr[j*4 + n+1];
			}
		}
	}
	for (int j = 0; j < nCells; j++) {
		vacc = 0;
		xi = xarr[j*4 +3];
		yi = yarr[j*4 +3];
		std::set<int> *nn = nnarr[j];
		for (std::set<int>::iterator it=nn->begin(); it!=nn->end(); ++it) {
			pn = *it;
			xn = xarr[pn*4 + 3];
			yn = yarr[pn*4 + 3];
			dx = xi - xn;
			dy = yi - yn;
			d = sqrt(dx*dx + dy*dy);
			vacc += 0.5*k*(d - d0)*(d - d0);
		}
		varr[j] = vacc;
	}
}

void Solver::initCells(double dt, double k, double l, double d0,
		std::vector<double>& xarr, std::vector<double>& yarr,
		std::vector<double>& vxarr, std::vector<double>& vyarr,
		std::vector<double>& fxarr, std::vector<double>& fyarr,
		std::vector<std::set<int>*>& nnarr, std::vector<double>& varr,
		std::vector<bool>& fiarr) {
	/*
	Create the initial values for the solver
	using euler method. The first value is left
	as the initial condition
	*/
	int nCells = fiarr.size();
	for (int j = 0; j < nCells; j++) {
		if (fiarr[j]) continue;
		calcAcc(0, 0, j, dt, k, l, d0, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
	}
	for (int i = 1; i < 4; i++) {
		for (int j = 0; j < nCells; j++) {
			if (fiarr[j]) continue;
			vxarr[j*4 + i] = vxarr[j*4 + i-1] + dt*fxarr[j*4 + i-1];
			vyarr[j*4 + i] = vyarr[j*4 + i-1] + dt*fxarr[j*4 + i-1];
			xarr[j*4 + i] = xarr[j*4 + i-1] + dt*vxarr[j*4 + i-1];
			yarr[j*4 + i] = yarr[j*4 + i-1] + dt*vyarr[j*4 + i-1];
			calcAcc(i, i-1, j, dt, k, l, d0, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
		}
	}
}

void Solver::initFreshCell(int numNewCells, double dt, double k, double l,
		double d0, std::vector<double>& xarr, std::vector<double>& yarr,
		std::vector<double>& vxarr, std::vector<double>& vyarr,
		std::vector<double>& fxarr, std::vector<double>& fyarr,
		std::vector<std::set<int>*>& nnarr, std::vector<double>& varr,
		std::vector<bool>& fiarr) {
	/*
	Create the previous 3 values for a new cell
	by going back in time
	*/

	for(unsigned int j = fiarr.size() - numNewCells; j < fiarr.size(); j++) {
		calcAcc(3, 2, j, dt, k, l, d0, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
		vxarr[j*4 + 2] = vxarr[j*4 + 2+1] - dt*fxarr[j*4 + 2+1];
		vyarr[j*4 + 2] = vyarr[j*4 + 2+1] - dt*fxarr[j*4 + 2+1];
		xarr[j*4 + 2] = xarr[j*4 + 2+1] - dt*vxarr[j*4 + 2+1];
		yarr[j*4 + 2] = yarr[j*4 + 2+1] - dt*vyarr[j*4 + 2+1];

		calcAcc(2, 1, j, dt, k, l, d0, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
		vxarr[j*4 + 1] = vxarr[j*4 + 1+1] - dt*fxarr[j*4 + 1+1];
		vyarr[j*4 + 1] = vyarr[j*4 + 1+1] - dt*fxarr[j*4 + 1+1];
		xarr[j*4 + 1] = xarr[j*4 + 1+1] - dt*vxarr[j*4 + 1+1];
		yarr[j*4 + 1] = yarr[j*4 + 1+1] - dt*vyarr[j*4 + 1+1];

		calcAcc(1, 0, j, dt, k, l, d0, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
		vxarr[j*4 + 0] = vxarr[j*4 + 1] - dt*fxarr[j*4 + 1];
		vyarr[j*4 + 0] = vyarr[j*4 + 1] - dt*fxarr[j*4 + 1];
		xarr[j*4 + 0] = xarr[j*4 + 1] - dt*vxarr[j*4 + 1];
		yarr[j*4 + 0] = yarr[j*4 + 1] - dt*vyarr[j*4 + 1];

		calcAcc(0, 0, j, dt, k, l, d0, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
	}
}

inline void Solver::calcAcc(int pa, int pb, int j, double dt, double k, double l,
		double d0, std::vector<double>& xarr, std::vector<double>& yarr,
		std::vector<double>& vxarr, std::vector<double>& vyarr,
		std::vector<double>& fxarr, std::vector<double>& fyarr,
		std::vector<std::set<int>*>& nnarr, std::vector<double>& varr,
		std::vector<bool>& fiarr) {
	double xi,yi,dx,dy,d, xn, yn;
	unsigned int pn;
	double acc[2];
	acc[0] = 0;
	acc[1] = 0;
	xi = xarr[j*4 + pa];
	yi = yarr[j*4 + pa];
	std::set<int> *nn = nnarr[j];
	for (std::set<int>::iterator it=nn->begin(); it!=nn->end(); ++it) {
		pn = *it;
		xn = xarr[pn*4 + pb];
		yn = yarr[pn*4 + pb];
		dx = xi - xn;
		dy = yi - yn;
		d = sqrt(dx*dx + dy*dy);
		acc[0] -= k*dx*(1-d0/d);
		acc[1] -= k*dy*(1-d0/d);
	}
	acc[0] -= l*vxarr[j*4 + pa];
	acc[1] -= l*vyarr[j*4 + pa];
	fxarr[j*4 + pa] = acc[0];
	fyarr[j*4 + pa] = acc[1];
}



