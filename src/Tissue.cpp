/*
 * Tissue.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: tiago
 */

#include "Tissue.h"

Tissue::Tissue(double k, double l, double d, int nx, int ny, int ns, int ruleID, double mx, double my) {
	_nx = nx;
	_ny = ny;
	_k = k;
	_l = l;
	_d = d;
	_pt = 0.5*k*pow(0.5*d,2);
	_md = 2*d;
	_hd = d/2;
	_eqTol = d*1e-2;
	createTransitionTables(ruleID);
	_dt = 0.01;
	_numNewCells = _nx*_ny;
	_numDead = 0;
	_avEquilib = 0;
	_avSize = 0;
	_numSteps = ns;
	_nIters = 0;
	allocateCells(mx, my);
	_mx = mx;
	_my = my*4;

	std::stringstream app;
	//TODO: make sure data folder exists, or else this crashes
	app << "/data/em" << getpid() << time(NULL) << ruleID << "/";
	_dataFolder = fs::current_path();
	_dataFolder /= app.str();
	std::cout << _dataFolder << std::endl;
	fs::create_directory(_dataFolder);
	saveCellState(0);

	std::cout << "Tissue setup for rule " << ruleID << std::endl;
}

Tissue::~Tissue() {
	for(std::vector<std::set<int>*>::iterator it = nnarr.begin(); it != nnarr.end(); it++)
		delete *it;
}

void Tissue::createTransitionTables(int ruleID) {
	transition.assign(8, 0);
	growth.assign(8, 0);
	_rid = ruleID;
	for (int i = 0; i < 8; i++) {
		transition[i] = (ruleID & 1);
		ruleID >>= 1;
	}
	for (int i = 0; i < 8; i++) {
		growth[i] = (ruleID & 1);
		ruleID >>= 1;
	}
	std::cout << "Rule " << _rid << " : ";
	for (int i = 0; i < 8; i++) std::cout << transition[i];
	for (int i = 0; i < 8; i++) std::cout << growth[i];
	std::cout << std::endl;
}

int Tissue::getRuleID() {
	return _rid;
}

void Tissue::evolveCA() {
	newCellX.clear();
	newCellY.clear();
	for (int j = 0; j < _nCells; j++) {
		int acc = 0;
		std::set<int> *nn = nnarr[j];
		for (std::set<int>::iterator it=nn->begin(); it!=nn->end(); ++it)
			acc += stateA[*it];
		int s = acc < 3 ? 0 : 1;
		int p = varr[j] < _pt ? 0 : 1;
		int x = stateA[j];
		int i = s<<2 | x<<1 | p;
		stateB[j] = transition[i];
		if (growth[i] == 1) {
			newCellX.push_back(xarr[j*4 + 3]);
			newCellY.push_back(yarr[j*4 + 3]);
		}
	}
	stateA.swap(stateB);
}

void Tissue::run() {
	clock_t begin = clock();
	for (int n = 0; n < _numSteps; n++) {
		/*
		Step 1: reach elastic equilibrium
		*/
		if (_numNewCells > 0 || _numDead > 0)
			equilibrate();

		/*
		Step 2: save current system state
		*/
		saveCellState(n);

		/*
		Step 3: update state and add cells
		*/
		evolveCA();
		_numNewCells = newCellX.size();

		/*
		Step 4: add cells (tCells)
		*/
		grow();
		_nCells = fiarr.size();

		/*
		Step 5: delete cells which went outside the tissue (tCells)
		*/
		_numDead = deleteRogueCells();
		_nCells = fiarr.size();
		if (_nCells > 1000) break;
	}
	equilibrate();
	clock_t end = clock();
	_rt = double(end - begin) / CLOCKS_PER_SEC;
	saveCellState(_numSteps);
	saveParameters();
}

void Tissue::equilibrate() {
	int i = 0;
	while(true) {
		i++;
		calcDelaunay();
		Solver::run(_dt, 2.0, _k, _l, _d, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
		double acc = 0;
		for(int i=0; i < _nCells; i++) {
			acc += sqrt((xarr[i*4 + 3] - xvert[i])*(xarr[i*4 + 3] - xvert[i]));
			acc += sqrt((yarr[i*4 + 3] - yvert[i])*(yarr[i*4 + 3] - yvert[i]));
		}
		acc /= 2.0 * _nCells;
		if (acc < _eqTol) break;
		else acc = 0;
	}
	std::cout << "equilib:" << i << " " << _nCells << std::endl;
	_avSize += _nCells;
	_avEquilib += i;
	_nIters += 1;


}

void Tissue::calcDelaunay() {
	xvert.resize(_nCells);
	yvert.resize(_nCells);
	for(int i = 0; i<_nCells; i++) {
		xvert[i] = xarr[i*4 + 3];
		yvert[i] = yarr[i*4 + 3];
		nnarr[i]->clear();
	}
	Triangulate::run(xvert, yvert, nnarr, _md);
}

void Tissue::grow() {
	for (int i = 0; i<_numNewCells; i++) {
		double rx = _hd*(mt_drand() - 0.5);
		double ry = _hd*(mt_drand() - 0.5);
		xarr.insert(xarr.end(), 4, newCellX[i] + rx);
		yarr.insert(yarr.end(), 4, newCellY[i] + ry);
		vxarr.insert(vxarr.end(), 4, 0);
		vyarr.insert(vyarr.end(), 4, 0);
		fxarr.insert(fxarr.end(), 4, 0);
		fyarr.insert(fyarr.end(), 4, 0);
		fiarr.push_back(false);
		varr.push_back(0);
		stateA.push_back(0);
		stateB.push_back(0);
		std::set<int> *nset = new std::set<int>;
		nnarr.push_back(nset);
	}
	Solver::initFreshCell(_numNewCells, _dt, _k, _l, _d, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
}

void Tissue::saveParameters() {
	fs::ofstream file(_dataFolder / "params.txt");
	//some c ugliness
	char buffer[512];
	time_t rawtime;
	time(&rawtime);
	struct tm * timeinfo = localtime(&rawtime);
	strftime(buffer, 512, "%d %b %Y %H:%M:%S", timeinfo);
	file << "Run finished at: " << buffer << std::endl;
	file << "Automaton rule: " << _rid << std::endl;
	file << "init num cells x: " << _nx << std::endl;
	file << "init num cells y: " << _ny << std::endl;
	file << "pressure threshold: " << _pt << std::endl;
	file << "lambda: "  << _l << std::endl;
	file << "stiffness: "  << _k << std::endl;
	file << "rest length: "  << _d << std::endl;
	file << "dt: " << _dt << std::endl;
	file << "number of ca steps: " << _numSteps << std::endl;
	file << "max x: " << _mx << std::endl;
	file << "max y: " << _my << std::endl;
	file << "average equilibration time: " << _avEquilib/(double)_nIters << std::endl;
	file << "average system size: " << _avSize/(double)_nIters << std::endl;
	file << "final system size: " << _nCells << std::endl;
	file << "runtime in seconds: " << _rt << std::endl;
	file << "state entropy: " << getStateEntropy() << std::endl;
	file << "network entropy: " << getNetworkEntropy() << std::endl;
	file.close();
}

void Tissue::printParameters() {
	std::cout << "Automaton rule: " << _rid << std::endl;
	std::cout << "init num cells x: " << _nx << std::endl;
	std::cout << "init num cells y: " << _ny << std::endl;
	std::cout << "pressure threshold: " << _pt << std::endl;
	std::cout << "lambda: "  << _l << std::endl;
	std::cout << "stiffness: "  << _k << std::endl;
	std::cout << "rest length: "  << _d << std::endl;
	std::cout << "dt: " << _dt << std::endl;
	std::cout << "number of ca steps: " << _numSteps << std::endl;
	std::cout << "max x: " << _mx << std::endl;
	std::cout << "max y: " << _my << std::endl;
	std::cout << "average equilibration time: " << _avEquilib/(double)_nIters << std::endl;
	std::cout << "average system size: " << _avSize/(double)_nIters << std::endl;
	std::cout << "final system size: " << _nCells << std::endl;
	std::cout << "runtime in seconds: " << _rt << std::endl;
	std::cout << "state entropy: " << getStateEntropy() << std::endl;
	std::cout << "network entropy: " << getNetworkEntropy() << std::endl;
}

void Tissue::saveCellState(int n) {
	char buffer[256];
	sprintf(buffer, "state%04d.cs", n);
	fs::ofstream file(_dataFolder / buffer);
	for(int i = 0; i < _nCells; i++) {
		file << xarr[i*4 + 3] << "; ";
	}
	file << std::endl;
	for(int i = 0; i < _nCells; i++) {
		file << yarr[i*4 + 3] << "; ";
	}
	file << std::endl;
	for(int i = 0; i < _nCells; i++) {
		file << stateA[i] << "; ";
	}
	file << std::endl;
	for(int i = 0; i < _nCells; i++) {
		file << varr[i] << "; ";
	}
		file << std::endl;
	for(int i = 0; i < _nCells; i++) {
		file << "[";
		std::set<int> *nn = nnarr[i];
		for (std::set<int>::iterator it=nn->begin(); it!=nn->end(); ++it) {
			file << *it << ",";
		}
		file << "]; ";
	}
	file << std::endl;
}

int Tissue::deleteRogueCells() {
	double tol = 1e-3;
	double bx = 0 - tol;
	double tx = 1 + tol;
	double by = 0 - tol;
	double ty = _my + tol;
	std::set<int> goners;

	for (int i = 0; i < _nCells; i++) {
		if (xarr[4*i + 3] < bx || xarr[4*i + 3] > tx
				|| yarr[4*i + 3] < by || yarr[4*i + 3] > ty) {
			goners.insert(i);
		}
	}

	for (std::set<int>::reverse_iterator it=goners.rbegin(); it!=goners.rend(); ++it) {
		int deadPos = *it;
		xarr.erase(xarr.begin() + 4*deadPos, xarr.begin() + 4*(deadPos + 1));
		yarr.erase(yarr.begin() + 4*deadPos, yarr.begin() + 4*(deadPos + 1));
		vxarr.erase(vxarr.begin() + 4*deadPos, vxarr.begin() + 4*(deadPos + 1));
		vyarr.erase(vyarr.begin() + 4*deadPos, vyarr.begin() + 4*(deadPos + 1));
		fxarr.erase(fxarr.begin() + 4*deadPos, fxarr.begin() + 4*(deadPos + 1));
		fyarr.erase(fyarr.begin() + 4*deadPos, fyarr.begin() + 4*(deadPos + 1));
		stateA.erase(stateA.begin() + deadPos);
		stateB.erase(stateB.begin() + deadPos);
		fiarr.erase(fiarr.begin() + deadPos);
		varr.erase(varr.begin() + deadPos);
		delete nnarr[deadPos];
		nnarr.erase(nnarr.begin() + deadPos);
	}

	return goners.size();
}

void Tissue::allocateCells(double mx, double my) {
	//make sure we arent always allocating
	xarr.reserve(4 * _numNewCells);
	yarr.reserve(4 * _numNewCells);
	fiarr.reserve(_numNewCells);
	for (int x = 0; x < _nx; x++) {
		for (int y = 0; y < _ny; y++) {
			double xPos = x/(double)(_nx-1.0)*mx + mt_drand()*1e-5;
			double yPos = y/(double)(_ny-1.0)*my + mt_drand()*1e-5;
			//or y == _ny-1;
			if(x == 0 || y == 0 || x == _nx-1)	{
				fiarr.push_back(true);
			}
			else {
				fiarr.push_back(false);
			}
			xarr.insert(xarr.end(), 4, xPos);
			yarr.insert(yarr.end(), 4, yPos);
		}
	}

	vxarr.assign(4 * _numNewCells, 0);
	vyarr.assign(4 * _numNewCells, 0);
	fxarr.assign(4 * _numNewCells, 0);
	fyarr.assign(4 * _numNewCells, 0);
	varr.assign(_numNewCells, 0);
	stateA.assign(_numNewCells, 0);
	stateB.assign(_numNewCells, 0);
	for(int i = 0; i < _numNewCells; i++) {
		std::set<int> *nset = new std::set<int>;
		nnarr.push_back(nset);
	}

	_nCells = _numNewCells;
	calcDelaunay();
	Solver::initCells(_dt, _k, _l, _d, xarr, yarr, vxarr, vyarr, fxarr, fyarr, nnarr, varr, fiarr);
}

double Tissue::getStateEntropy() {
	double acc;
	for(int i=0; i<_nCells; i++) {
		acc += stateA[i];
	}
	double f1 = acc/(double)_nCells;
	double f0 = 1 - f1;
	return xlogx(f1) + xlogx(f0);
}


double Tissue::getNetworkEntropy() {
	std::vector<double> lenFreqs(10, 0);
	for(int i=0; i<_nCells; i++) {
		lenFreqs[nnarr[i]->size()] += 1;
	}
	double acc = 0;
	for(std::vector<double>::iterator it=lenFreqs.begin(); it!=lenFreqs.end(); ++it) {
		double fi = *it/(double)_nCells;
		acc += xlogx(fi);
	}
	return acc;
}

inline double Tissue::xlogx(double x) {
	double lx = x > 0 ? log(x) : 0;
	return x*lx;
}


