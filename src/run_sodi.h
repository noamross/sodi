
#ifndef RUNSODI
#define RUNSODI

#include "includes.h"
#include "data_structures.h"

// extern int run_sodi_rcpp(DataFrame init, List parm, bool progress, CharacterVector file, bool diagnostics, CharacterVector diagname);
extern double (*beta_f)(double beta, int i, int max_inf);
extern arma::vec (*kernel1)(arma::vec d, arma::vec m);
extern arma::vec (*kernel2)(arma::vec d, double m);
extern arma::vec (*seedkern)(double x, double y, int s, parmlist &parms);
extern void mgmt_action(statelist &state, parmlist &parms, arma::mat R);

#endif /* RUNSODI */