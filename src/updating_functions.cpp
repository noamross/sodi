#include "includes.h"
#include "data_structures.h"
#include "run_sodi.h"
#include "dispersal_functions.h"
#include "lamda_interp.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

void reproduce(statelist &state, parmlist &parms, arma::mat &R, uword j, uword s, uword treecount) {
  state.IDcount += 1;
  state.ID(treecount) = state.IDcount;
  arma::vec loc = seedkern(state.X(j), state.Y(j), s, parms);
  state.X(treecount) = loc(0);
  state.Y(treecount) = loc(1);
  state.S(treecount) = parms.ss(s);
  state.I(treecount) = 0;
  state.B(treecount) = parms.beta(s);
  state.F(treecount) = sum(kernel1(distance(state, treecount), parms.m(state.S)) % state.I % parms.lamda(state.S));
  state.E = fmax(0, 1 - sum(parms.omega(state.S)) / parms.K);
  R.col(0) = state.E * parms.f(state.S);
  R(treecount, 1) = parms.d(parms.ss(s));
  R(treecount, 2) = parms.g(parms.ss(s));
  R(treecount, 3) = (state.F(treecount) + lamda_interp(state, parms)) * state.B(treecount);
  R(treecount, 4) = 0;
  R(treecount, 5) = 0;
  state.treecount += 1;
  state.treeindex += 1;
}

void die(statelist &state, parmlist &parms, arma::mat &R, uword j, uword s, uword treeindex) {
   state.F -= kernel2(distance(state, j), parms.m(s))  * state.I(j) * parms.lamda(s);
   for (arma::uword z = 0; z < state.treecount; z++) {
      state.F(z) = fmax(0, state.F(z));
    }
   state.ID(j) = state.ID(treeindex);
   state.X(j) = state.X(treeindex);
   state.Y(j) = state.Y(treeindex);
   state.S(j) = state.S(treeindex);
   state.I(j) = state.I(treeindex);
   state.B(j) = state.B(treeindex);
   state.F(j) = state.F(treeindex);
   R.row(j) = R.row(treeindex);
   state.ID(treeindex) = 0;
   state.X(treeindex) = 0;
   state.Y(treeindex) = 0;
   state.S(treeindex) = 0;
   state.I(treeindex) = 0;
   state.B(treeindex) = 0;
   state.F(treeindex) = 0;
   R.row(treeindex).zeros();
   state.E = fmax(0, 1 - sum(parms.omega(state.S)) / parms.K);
   R.col(0) = state.E * parms.f(state.S);
   R.col(3) = (state.F + lamda_interp(state, parms)) % state.B;    
   state.treecount -= 1;
   state.treeindex -= 1;
}

void grow(statelist &state, parmlist &parms, arma::mat &R, uword j, uword s, int i) {

   state.S(j) += 1;
   state.B(j) = beta_f(parms.beta(state.S(j)), i, parms.max_inf(state.S(j)));
   state.F += kernel2(distance(state, j), parms.m(state.S(j))) * i * (parms.lamda(state.S(j)) - parms.lamda(s));
   state.E = fmax(0, 1 - sum(parms.omega(state.S)) / parms.K);
   R.col(0) = state.E * parms.f(state.S);
   R(j, 1) = parms.d(state.S(j)) + i * parms.alpha(state.S(j)) * (1 - parms.r(state.S(j)));
   R(j, 2) = parms.g(state.S(j));
   R.col(3) = (state.F + lamda_interp(state, parms)) % state.B;
   R(j, 4) = i * parms.mu(state.S(j));
   R(j, 5) = i * parms.alpha(state.S(j)) * parms.r(state.S(j));
}
    
void sicken(statelist &state, parmlist &parms, arma::mat &R, uword j, uword s, int i) {
   state.I(j) += 1;
   state.F += kernel2(distance(state, j), parms.m(s)) * parms.lamda(s);
   state.B(j) = beta_f(parms.beta(s), state.I(j), parms.max_inf(s));
   state.F(j) = state.F(j) * state.B(j) / beta_f(parms.beta(s), i, parms.max_inf(s));
   R(j, 1) += parms.alpha(s) * (1 - parms.r(s));
   R.col(3) = (state.F + lamda_interp(state, parms)) % state.B;
   R(j, 4) += parms.mu(s);
   R(j, 5) += parms.alpha(s) * parms.r(s);
   }

void recover(statelist &state, parmlist &parms, arma::mat &R, uword j, uword s, int i) {
   state.I(j) -= 1;
   state.F -= kernel2(distance(state, j), parms.m(s))  * parms.lamda(s);
   for (arma::uword z = 0; z < state.treecount; z++) {
     state.F(z) = fmax(0, state.F(z));
   }
   state.B(j) = beta_f(parms.beta(s), state.I(j), parms.max_inf(s));
   R(j, 1) -= parms.alpha(s) * (1 - parms.r(s));
   R(j, 1) = fmax(0, R(j, 1));
   R.col(3) = (state.F + lamda_interp(state, parms)) % state.B;
   R(j, 4) -= parms.mu(s);
   R(j, 4) = fmax(0, R(j, 4));
   R(j, 5) -= parms.alpha(s) * parms.r(s);
   R(j, 5) = fmax(0, R(j, 5));
   }

void resprout(statelist &state, parmlist &parms, arma::mat &R, uword j, uword s, int i) {
   state.S(j) = parms.ss(s);  //Set the stage back to the "base stage" for the species
   state.I(j) = 0;
   state.B(j) = parms.beta((state.S(j)));
   state.F -= kernel2(distance(state, j), parms.m(s)) * i * parms.lamda(s);
   state.F(j) = sum(kernel1(distance(state, j), parms.m(state.S)) % state.I % parms.lamda(state.S)) * state.B(j);
   for (arma::uword z = 0; z < state.treecount; z++) {
     state.F(z) = fmax(0, state.F(z));
   }
   R(j, 0) = parms.f(state.S(j));
   R(j, 1) = parms.d(state.S(j));
   R(j, 2) = parms.g(state.S(j));
   R.col(3) = (state.F + lamda_interp(state, parms)) % state.B;
   R(j, 4) = 0;
   R(j, 5) = 0;
}