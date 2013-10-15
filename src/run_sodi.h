
#ifndef RUNSODI
#define RUNSODI

#include <RcppArmadillo.h>
#include <fstream>
#include <string>
#include "data_structures.h"
#include "dispersal_functions.h"
#include "infection_density_dependence.h"
#include "updating_functions.h"
#include "print_state.h"
#include "data_structures.h"

extern int run_sodi_rcpp(Rcpp::DataFrame init, Rcpp::List parm, bool progress, Rcpp::CharacterVector file);
extern double (*beta_f)(double beta, int i, int max_inf);
extern arma::vec (*kernel1)(arma::vec d, arma::vec m);
extern arma::vec (*kernel2)(arma::vec d, double m);
extern arma::vec (*seedkern)(double x, double y, int s, parmlist &parms);

#endif /* RUNSODI */