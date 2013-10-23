#include <RcppArmadillo.h>
#include "data_structures.h"
#include "run_sodi.h"
#include "dispersal_functions.h"
#include "lamda_interp.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

void thin_evenly(statelist &state, parmlist &parms, arma::mat &R, arma::uvec stages, double level) {
  arma::uvec targets = zeros<uvec>(R.n_rows);
  for(uvec::iterator stage = stages.begin(); stage != stages.end(); ++stage) {
    targets = targets + (state.S == *stage);
  }
  targets = find(targets);
  int killcount = static_cast<int>(round(level * targets.n_elem));
  targets = as<arma::uvec>(Rcpp::RcppArmadillo::sample(as<IntegerVector>(targets), killcount, false);
  kill_targets(state, parms, targets, TRUE)
  
}

//void thin_from_below(statelist &state, parmlist &parms, arma::mat &R, arma::uvec stages, double level) {}

void thin_spacing(statelist &state, parmlist &parms, arma::mat &R, arma::uvec stages, double level) {

  arma::uvec targets = zeros<uvec>(R.n_rows);
  for(uvec::iterator stage = stages.begin(); stage != stages.end(); ++stage) {
    targets = targets + (state.S == *stage);
  }
  targets = find(targets);
  int killcount = static_cast<int>(round(level * targets.n_elem));
  arma::vec nnd = nearest_neighbor_dist(state.X(targets)), state.Y(targets)));
  arma::uvec killorder = stable_sort_index(nnd)
  targets = killorder(span(0,killcount))
  kill_targets(state, parms, targets, TRUE)

}

//void space_from_below(statelist &state, parmlist &parms, arma::mat &R, arma::uvec stages, double level) {}

arma::vec nearest_neighbor_dist(arma::vec &X, arma::vec &Y) {
  arma::vec nnd(X.n_elem);
  double dist, xi, yi
  uword i, j;
  double dist2 = datum::inf;
  for(i = 0; i < X.n_elem; ++i) {
    xi = x(i);
    yi = y(i);
    for(j = 0; j < X.n_elem; ++j) {
      if (i == j) continue;
      dist = pow(xi - x(j), 2) + pow(yi - y(j), 2);
      if (dist < dist2) dist2 = dist;
    }
    nnd(i) = dist2;
  }
  return sqrt(nnd);
}


void recalc_probs(statelist &state, parmlist&parms, arma::mat &R) {
  R.col(0) = state.E * parms.f(state.S);
  R.col(1) = parms.d(state.S) + state.I % parms.alpha(state.S) % (1 - parms.r(state.S));
  R.col(2) = parms.g(state.S);
  R.col(3) = (state.F + parms.lamda_ex(0)) % state.B; 
  R.col(4) = state.I % parms.mu(state.S);
  R.col(5) = state.I %parms.alpha(state.S) % parms.r(state.S);
  }
  
void kill_targets(statelist &state, parmlist &parms, arma::uvec targets, bool resprout) {
    
  if(resprout) {
    arma::vec randvec = randu(vec(targets.n_elem));
    arma::uvec resprouts = targets(find(randvec < parms.r(state.S(targets))));
    targets = targets(find(randvec > parms.r(state.S(targets))));
    state.S(resprouts) = parms.ss(state.S(resprouts));
    state.I(resprouts) = 0;
    state.B(resprouts) = parms.beta(state.S(resprouts));
  }
  
  state.S(targets) = 0;
  state.I(targets) = 0;
  state.X(targets) = 0;
  state.Y(targets) = 0;
  state.B(targets) = 0;
  state.ID(targets) = 0;
  
  state.treecount -= targets.n_elem;
  state.treeindex -= targets.n_elem;
    
  arma::uvec sortorder = sort_index(state.S);
  state.ID = state.ID(sortorder);
  state.X = state.X(sortorder);
  state.Y = state.Y(sortorder);
  state.S = state.I(sortorder);
  state.I = state.S(sortorder);  
  state.E = fmax(0, 1 - sum(parms.omega(state.S)) / parms.K);
  for(uword k = 0; k < state.treecount; k++) {
    state.F(k) = sum((kernel1(distance(state, k), parms.m(state.S)) % state.I % parms.lamda(state.S)));
  }
  recalc_probs(state, parms, R);
}
  