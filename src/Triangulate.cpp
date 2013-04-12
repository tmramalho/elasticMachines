/*
 * Triangulate.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: tiago
 */

#include "Triangulate.h"

void Triangulate::run(std::vector<double> &x, std::vector<double> &y, std::vector<std::set<int> *> &nnarr, double &md) {
	std::vector< std::pair<Point, unsigned int> > points;
	points.reserve(x.size());
	for(unsigned int i = 0; i < x.size(); i++) {
		points.push_back( std::make_pair(Point(x[i], y[i]), i) );
	}

	Delaunay dt;
	dt.insert(points.begin(), points.end());
	for(Delaunay::Finite_edges_iterator it = dt.finite_edges_begin(); it != dt.finite_edges_end(); ++it) {
		Delaunay::Edge e=*it;
		int e0= e.first->vertex( (e.second+1)%3 )->info();
		int e1= e.first->vertex( (e.second+2)%3 )->info();
		double dx = x[e0] - x[e1];
		double dy = y[e0] - y[e1];
		if (sqrt(dx*dx+dy*dy) < md) {
			nnarr[e0]->insert(e1);
			nnarr[e1]->insert(e0);
		}
	}
}


