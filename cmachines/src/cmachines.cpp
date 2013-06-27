//============================================================================
// Name        : C Machines
// Author      : Tiago
// Version     : 1
// Copyright   :
// Description :
//============================================================================

#include "mtwist.h"
#include <gflags/gflags.h>

#include <iostream>
#include <time.h>
#include <vector>
#include <string>
#include <unistd.h>
#include "Tissue.h"
#include "omp.h"
#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

DEFINE_double(k, 3.0, "elastic constant for the cells");
DEFINE_double(l, 1.0, "dampening constant");
DEFINE_double(d, 0.01, "equilibrium intra nuclei distance");
DEFINE_int32(nx, 100, "initial number of cells in x dim");
DEFINE_int32(ny, 3, "initial number of cells in x dim");
DEFINE_int32(ns, 10, "number of cellular automata iterations");
DEFINE_int32(rid, 19574, "integer describing the cellular automata transition table.Must be < 2^16, else mod will be taken");
DEFINE_int32(scb, 0, "Scan a range of rules[Begin], if selected --rid will be ignored.");
DEFINE_int32(sce, 0, "Scan a range of rules[End], if selected --rid will be ignored.");
DEFINE_int32(nth, 1, "Number of parallel threads. Only makes sense if scanning over a range of rules.");
DEFINE_double(xe, 1.0, "Cells will be initially uniformely distributed with positions in range [0, x_end]");
DEFINE_double(ye, 0.02, "Cells will be initially uniformely distributed with positions in range [0, x_end]");

int main(int argc, char *argv[]) {

	google::ParseCommandLineFlags(&argc, &argv, true);
	mt_seed32( time(NULL) + getpid() );

	fs::path baseDataFolder = fs::current_path() / "data";
	if(!fs::is_directory(baseDataFolder)) {
		if (exists(baseDataFolder)) {
			std::cerr << "I need to be able to access the /data folder." << std::endl;
			exit(1);
		} else {
			fs::create_directory(baseDataFolder);
		}
	}

	if(FLAGS_scb == 0 && FLAGS_sce == 0) {
		Tissue t(FLAGS_k, FLAGS_l, FLAGS_d, FLAGS_nx, FLAGS_ny, FLAGS_ns, FLAGS_rid, FLAGS_xe, FLAGS_ye);
		t.run();
		t.writeParameters(std::cout);
	} else {
		omp_set_num_threads(FLAGS_nth);
		if(FLAGS_sce <= FLAGS_scb || FLAGS_scb < 0) {
			std::cerr << "Range spec wrong: " << FLAGS_scb << " to " << FLAGS_sce << std::endl;
			exit(1);
		}
		#pragma omp parallel for
		for(int i = FLAGS_scb; i <= FLAGS_sce; i++) {
			Tissue t(FLAGS_k, FLAGS_l, FLAGS_d, FLAGS_nx, FLAGS_ny, FLAGS_ns, i, FLAGS_xe, FLAGS_ye);
			t.run();
		}
	}

	return 0;
}
