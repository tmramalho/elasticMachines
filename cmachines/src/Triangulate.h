/*
 * Triangulate.h
 *
 *  Created on: Apr 10, 2013
 *      Author: tiago
 */

#ifndef TRIANGULATE_H_
#define TRIANGULATE_H_
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                   		Delaunay;
typedef K::Point_2                                          		Point;

class Triangulate {
public:
	static void run(std::vector<double> &x, std::vector<double> &y, std::vector<std::set<int> *> &nnarr, double &md);
};

#endif /* TRIANGULATE_H_ */
