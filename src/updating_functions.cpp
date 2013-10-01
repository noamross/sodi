#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
void reproduce(statelist& state, parmlist& parms, uword j, uword s, uword treecount) {
  state.IDcount += 1;
  state.ID(treecount) = state.IDcount;
  arma::vec loc = seedkern(X(j), Y(j), s, &parms);
  state.X(treecount) = loc(1);
  state.Y(treecount) = loc(2);
  state.S(treecount) = parms.ss(s);
  state.I(treecount) = 0;
  state.B(treecount) = parms.beta(s);
  state.F(treecount) = sum(kernel(distance(&state, treecount), m(S)) % I % lamda(S)) * B(treecount);
  state.E = fmax(0, 1 - sum(parms.omega(state.S)) / parms.K);
  R.col(1) = state.E * parms.f(state.S);
  R(treecount, 2) = parms.d(parms.ss(s))
  R(treecount, 3) = parms.g(parms.ss(s));
  R(treecount, 4) = state.F(treecount);
  R(treecount, 5) = 0;
  R(treecount, 6) = 0;
  state.treecount += 1;
  state.treeindex += 1;
}

void die(statelist& state, parmlist& parms, uword j, uword s, uword treeindex) {
   state.F -= kernel2(distance(&state, j), m(s)) % state.B * i * lamda(s);
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
   R.col(1) = state.E * parms.f(state.S);
   R.col(4) = state.F; 
   state.treecount -= 1;
   state.treeindex -= 1;
}

void grow(statelist &state, parmlist &parms, uword j, uword s, int i) {

   state.S(j) += 1;
   state.B(j) = beta_f(parms.beta(state.S(j)), i, parms.max_inf(S(j)));
   state.F(j) = F(j) * parms.beta(state.S(j))/parms.beta(s);
   F += kernel2(distance(X, Y, j), m(S(j))) % B * i * (lamda(S(j)) - lamda(s));
   
   state.E = fmax(0, 1 - sum(parms.omega(state.S)) / parms.K);
   R.col(1) = state.E * parms.f(state.S);
   R(j, 2) = parms.d(state.S(j)) + i * alpha(state.S(j)) % (1 - parms.r(state.S(j)));
   R(j, 3) = parms.g(state.S(j));
   R.col(4) = state.F;
   R(j, 5) = i % parms.mu(state.S(j));
   R(j, 6) = i % alpha(state.S(j)) * parms.r(state.S(j));
}
    
void sicken(statelist &state, parmlist &parms, uword j, uword s, int i) {
   state.I(j) += 1;
   F += kernel2(distance(&state, j), parms.m(s)) % state.B * parms.lamda(s);
   state.B(j) = parms.beta_f(beta(s), state.I(j), parms.max_inf(s));
   F(j) = F(j) * B(j) / beta_f(beta(s), i, max_inf(s));
   R(j, 2) += parms.alpha(s) * (1 - parms.r(s));
   R.col(4) = state.F;
   R(j, 5) += parms.mu(s);
   R(j, 6) += parms.alpha(s) * parms.r(s);
   }

void recover(statelist &state, parmlist &parms, uword j, uword s, int i) {
   state.I(j) -= 1;
   F -= kernel2(distance(&state, j), parms.m(s)) % state.B * parms.lamda(s);
   state.B(j) = parms.beta_f(beta(s), state.I(j), parms.max_inf(s));
   F(j) = F(j) * B(j) / beta_f(beta(s), i, max_inf(s));
   R(j, 2) -= parms.alpha(s) * (1 - parms.r(s));
   R.col(4) = state.F;
   R(j, 5) -= parms.mu(s);
   R(j, 6) -= parms.alpha(s) * parms.r(s);
   }

void resprout(statelist &state, parmlist &parms, uword j, uword s, int i) {
   state.S(j) = parms.ss(s);  //Set the stage back to the "base stage" for the species
   I(j) = 0;
   B(j) = beta((S(j)));
   F -= kernel2(distance(&state, j), parms.m(s)) % state.B * i * parms.lamda(s);
   F(j) = sum(kernel(distance(&state, j), parms.m(S)) % state.I % parms.lamda(S)) * state.B(j);
   R(j, 1) = parms.f(state.S(j));
   R(j, 2) = parms.d(state.S(j));
   R(j, 3) = parms.g(state.S(j));
   R.col(4) = state.F;
   R(j, 5) = 0
   R(j, 6) = 0
}