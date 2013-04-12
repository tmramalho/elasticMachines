//============================================================================
// Name        : C Machines
// Author      : Tiago
// Version     : 1
// Copyright   :
// Description :
//============================================================================

#include "mtwist.h"

#include <iostream>
#include <time.h>
#include <vector>
#include <string>
#include <unistd.h>
#include "Tissue.h"

int main(int argc, char *argv[]) {

	mt_seed32( time(NULL) + getpid() );

	Tissue t(3.0, 1.0, 0.01, 100, 3, 10, 19574, 1.0, 0.02);
	t.run();
	t.printParameters();

	return 0;
}
