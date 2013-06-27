/*
 * Tissue.h
 *
 *  Created on: Apr 10, 2013
 *      Author: tiago
 */

#ifndef TISSUE_H_
#define TISSUE_H_
#include <iosfwd>
#include <cmath>
#include <iostream>
#include <vector>
#include <set>
#include <sstream>
#include <ctime>
#include <cmath>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
namespace fs = boost::filesystem;
#include "Solver.h"
#include "mtwist.h"
#include "Triangulate.h"

class Tissue {
public:
	Tissue(double k, double l, double d, int nx, int ny, int ns, int ruleID, double mx, double my);
	~Tissue();
	void createTransitionTables(int ruleID);
	int getRuleID();
	void run();
	void writeParameters(std::ostream& stream);
	double getStateEntropy();
	double getNetworkEntropy();

private:
	void equilibrate();
	void calcDelaunay();
	void grow();
	void evolveCA();
	int deleteRogueCells();
	void saveCellState(int n, std::ostream& stream);
	void allocateCells(double mx, double my);
	inline double xlogx(double x);
	int _nx, _ny, _numNewCells, _numDead, _numMaxSteps, _nIters, _rid, _nCells, _numTotalSteps;
	double _pt, _md, _hd, _eqTol, _dt, _avEquilib, _avSize, _mx, _my, _k, _l, _d, _rt;
	std::vector<double> xarr;
	std::vector<double> yarr;
	std::vector<double> xvert;
	std::vector<double> yvert;
	std::vector<double> vxarr;
	std::vector<double> vyarr;
	std::vector<double> fxarr;
	std::vector<double> fyarr;
	std::vector<std::set<int>*> nnarr;
	std::vector<double> varr;
	std::vector<bool> fiarr;
	std::vector<int> stateA;
	std::vector<int> stateB;
	std::vector<int> transition;
	std::vector<int> growth;
	std::vector<double> newCellX;
	std::vector<double> newCellY;
	fs::path _dataFolder;
};

#endif /* TISSUE_H_ */
