#ifndef PRINT
#define PRINT

#include "includes.h"
#include <fstream>
#include "data_structures.h"


extern void print_state(statelist &state, std::ofstream &outfile, arma::mat &printmatrix);

#endif /* PRINT */

