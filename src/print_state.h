#ifndef PRINT
#define PRINT
#include <RcppArmadillo.h>
#include "data_structures.h"

extern void print_state(double outtime, statelist &state, std::ofstream &outfile, arma::mat &printmatrix);
#endif