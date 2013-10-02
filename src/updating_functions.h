#ifndef UPFUN
#define UPFUN

#include <RcppArmadillo.h>
#include "data_structures.h"
#include "run_sodi.h"

extern void reproduce(statelist &state, parmlist &parms, arma::mat &R, arma::uword j, arma::uword s, arma::uword treecount);
extern void die(statelist &state, parmlist &parms, arma::mat &R, arma::uword j, arma::uword s, arma::uword treeindex);
extern void grow(statelist &state, parmlist &parms, arma::mat &R, arma::uword j, arma::uword s, int i);
extern void sicken(statelist &state, parmlist &parms, arma::mat &R, arma::uword j, arma::uword s, int i);
extern void recover(statelist &state, parmlist &parms, arma::mat &R, arma::uword j, arma::uword s, int i); 
extern void resprout(statelist &state, parmlist &parms, arma::mat &R, arma::uword j, arma::uword s, int i);

#endif