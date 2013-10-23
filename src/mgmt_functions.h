
#ifndef MGMT
#define MGMT
#include "includes.h"
#include "data_structures.h"

// extern int run_sodi_rcpp(DataFrame init, List parm, bool progress, CharacterVector file, bool diagnostics, CharacterVector diagname);
void thin_evenly(statelist &state, parmlist &parms, arma::mat &R, arma::uvec stages, double level, bool resprout);
void thin_spacing(statelist &state, parmlist &parms, arma::mat &R, arma::uvec stages, double level, bool resprout);
arma::vec nearest_neighbor_dist(arma::vec &X, arma::vec &Y);
void recalc_probs(statelist &state, parmlist&parms, arma::mat &R);
void kill_targets(statelist &state, parmlist &parms, arma::mat &R, arma::uvec targets, bool resprout);

#endif /* MGMT */