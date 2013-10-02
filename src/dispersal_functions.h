#ifndef DISPFUN
#define DISPFUN

#include <RcppArmadillo.h>
#include "data_structures.h"

extern arma::vec flatdisp(arma::vec d, arma::vec m);
extern arma::vec flatdisp2(arma::vec d, double m);
extern arma::vec expdisp(arma::vec d, arma::vec m);
extern arma::vec expdisp2(arma::vec d, double m);
extern arma::vec normdisp(arma::vec d, arma::vec m);
extern arma::vec normdisp2(arma::vec d, double m);
extern arma::vec fatdisp(arma::vec d, arma::vec m);
extern arma::vec fatdisp2(arma::vec d, double m);
extern arma::vec distance(statelist &state, int j);
extern arma::vec powvec(arma::vec A, arma::ivec p);
extern arma::vec flatseed(double x, double y, int s, parmlist &parms);
extern arma::vec normseed(double x, double y, int s, parmlist &parms);
extern arma::vec expseed(double x, double y, int s, parmlist &parms);



#endif