/*
 * Solver.h
 *
 *  Created on: Apr 10, 2013
 *      Author: tiago
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <vector>
#include <cmath>
#include <set>

class Solver {
public:
	static void run(double dt, double total, double k, double l, double d0,
					std::vector<double> &xarr,
					std::vector<double> &yarr,
					std::vector<double> &vxarr,
					std::vector<double> &vyarr,
					std::vector<double> &fxarr,
					std::vector<double> &fyarr,
					std::vector<std::set<int> *> &nnarr,
					std::vector<double> &varr,
					std::vector<bool> &fiarr);
	static void initCells(double dt, double k, double l, double d0,
					std::vector<double> &xarr,
					std::vector<double> &yarr,
					std::vector<double> &vxarr,
					std::vector<double> &vyarr,
					std::vector<double> &fxarr,
					std::vector<double> &fyarr,
					std::vector<std::set<int> *> &nnarr,
					std::vector<double> &varr,
					std::vector<bool> &fiarr);
	static void initFreshCell(int numNewCells, double dt, double k, double l, double d0,
					std::vector<double> &xarr,
					std::vector<double> &yarr,
					std::vector<double> &vxarr,
					std::vector<double> &vyarr,
					std::vector<double> &fxarr,
					std::vector<double> &fyarr,
					std::vector<std::set<int> *> &nnarr,
					std::vector<double> &varr,
					std::vector<bool> &fiarr);
	inline static void calcAcc(int pa, int pb, int j, double dt, double k, double l, double d0,
			std::vector<double> &xarr,
			std::vector<double> &yarr,
			std::vector<double> &vxarr,
			std::vector<double> &vyarr,
			std::vector<double> &fxarr,
			std::vector<double> &fyarr,
			std::vector<std::set<int> *> &nnarr,
			std::vector<double> &varr,
			std::vector<bool> &fiarr);
};

#endif /* SOLVER_H_ */
