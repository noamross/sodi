#ifndef STRUCTURES
#define STRUCTURES

#include <RcppArmadillo.h>

struct parmlist {
    int K;
    int N;
    arma::vec bbox;
    arma::vec seedm;
    arma::vec m;
    arma::vec f;
    arma::vec g;
    arma::vec d;
    arma::vec r;
    arma::vec alpha;
    arma::vec lamda;
    arma::vec beta;
    arma::vec mu;
    arma::vec xi;
    arma::vec omega;
    arma::uvec ss;
    arma::ivec max_inf;
    
  public:
};

struct statelist {
    arma::uvec ID;
    arma::vec X; 
  	arma::vec Y; 
    arma::uvec S;
    arma::ivec I;
    arma::vec F;
    arma::vec B;
    arma::vec R;
    double time;
    arma::uword IDcount;
    arma::uword treecount;  
    arma::uword treeindex;
    double E;
};


#endif /* STRUCTURES */