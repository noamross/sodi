#ifndef STRUCTURES
#define STRUCTURES

#include "includes.h"

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
    arma::vec lamda_ex;
    arma::vec times;
    arma::uvec mg_actions;
    arma::uvec mg_actionlist;
    arma::vec mg_levels;
    Rcpp::LogicalVector mg_resprout;
    arma::umat mg_stages;
    
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
    arma::vec::iterator next_record;
    arma::uvec::iterator next_act;
};


#endif /* STRUCTURES */