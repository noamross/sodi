#include <RcppArmadillo.h>
#include "data_structures.h"

//The dispersal functions - generate probability of infection from distance 
arma::vec flatdisp(arma::vec d, arma::vec m) {
  return m;
}

arma::vec flatdisp2(arma::vec d, double m) {
  d.fill(m);
  return d;
}

arma::vec expdisp(arma::vec d, arma::vec m) {
  return (square(m) / (2*M_PI)) % exp(-m % d);
}

arma::vec expdisp2(arma::vec d, double m) {
  return (m*m / (2*M_PI)) * exp(-m * d);
}

arma::vec normdisp(arma::vec d, arma::vec m) {
  return (square(m) / sqrt(2*M_PI)) % exp(-0.5 * square(m % d));
}

arma::vec normdisp2(arma::vec d, double m) {
  return (m*m / sqrt(2*M_PI)) * exp(-0.5 * square(m * d));
}

arma::vec fatdisp(arma::vec d, arma::vec m) {
  return (square(m) / (24* M_PI)) % exp(- sqrt(m % d));
}

arma::vec fatdisp2(arma::vec d, double m) {
  return (m*m / (24* M_PI)) * exp(- sqrt(m * d));
}


//Distance calculation.  Returns vector of Euclidian distances between tree j
//and the rest of trees
arma::vec distance(statelist &state, int j) {
  return sqrt(square(state.X - state.X(j)) + square(state.Y - state.Y(j)));
}

//Vectorized power calculation
arma::vec powvec(arma::vec A, arma::ivec p) {
  arma::uword n = A.n_elem;
  arma::vec B(n);
  for(arma::uword i = 0; i < n; ++i) {
    B(i) = pow(A(i), p(i));
  }
  return B;
}

//Seed disperal functions - generate random locations for new trees

arma::vec flatseed(double x, double y, int s, parmlist &parms) {
  arma::vec loc(2);
  loc(1) = Rcpp::as<double>(Rcpp::runif(1, parms.bbox(0), parms.bbox(1)));
  loc(2) = Rcpp::as<double>(Rcpp::runif(1, parms.bbox(2), parms.bbox(3)));
  return loc;
}

arma::vec normseed(double x, double y, int s, parmlist &parms) {
  arma::vec loc(2);
  double theta = Rcpp::as<double>(Rcpp::runif(1, 0, 2*M_PI));
  double dist = Rcpp::as<double>(Rcpp::rnorm(1, parms.seedm(s)));
  loc(1) = x + cos(theta) * dist;
  loc(2) = y + sin(theta) * dist;
  return loc;
}

arma::vec expseed(double x, double y, int s, parmlist &parms) {
  arma::vec loc(2);
  double theta = Rcpp::as<double>(Rcpp::runif(1, 0, 2*M_PI));
  double dist = Rcpp::as<double>(Rcpp::rexp(1, parms.seedm(s)));
  loc(1) = x + cos(theta) * dist;
  loc(2) = y + sin(theta) * dist;
  return loc;
}
