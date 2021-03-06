#define PROFILE 0
#define PROFILE_FILE "~/out.prof"

#if PROFILE
#include <gperftools/profiler.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include "includes.h"
#include "mgmt_functions.h"
#include <fstream>
#include <string>
#include "data_structures.h"
#include "dispersal_functions.h"
#include "infection_density_dependence.h"
#include "updating_functions.h"
#include "print_state.h"
#include "run_sodi.h"
#include "lamda_interp.h"


using namespace Rcpp;
using namespace arma;
using namespace std;


arma::vec (*kernel1)(arma::vec d, arma::vec m);
arma::vec (*kernel2)(arma::vec d, double m);
arma::vec (*seedkern)(double x, double y, int s, parmlist &parms);
double (*beta_f)(double beta, int i, int max_inf);


//' @export
// [[Rcpp::export]]
int run_sodi_rcpp(DataFrame init, List parm, bool progress, CharacterVector file, bool diagnostics, CharacterVector diagname) {
  
  #if PROFILE
  ProfilerStart(PROFILE_FILE);
  #endif
  
  std::string filename = Rcpp::as<std::string>(file);
  std::string diagfile = Rcpp::as<std::string>(diagname);
  
  //populate the parameter structure from the input list
  parmlist parms;
  parms.K = as<int>(parm["K"]);
  parms.N = as<int>(parm["n0"]);
  parms.bbox = as<NumericVector>(parm["bbox"]);
  parms.seedm = as<arma::vec>(parm["seedm"]);
  parms.m = as<arma::vec>(parm["m"]);
  parms.f = as<arma::vec>(parm["f"]);
  parms.g = as<arma::vec>(parm["g"]);
  parms.d = as<arma::vec>(parm["d"]);
  parms.r = as<arma::vec>(parm["r"]);
  parms.alpha = as<arma::vec>(parm["alpha"]);
  parms.lamda = as<arma::vec>(parm["lamda"]);
  parms.beta = as<arma::vec>(parm["beta"]);
  parms.mu = as<arma::vec>(parm["mu"]);
  parms.xi = as<arma::vec>(parm["xi"]);
  parms.omega = as<arma::vec>(parm["omega"]);
  parms.ss = as<arma::uvec>(parm["ss"]);
  parms.max_inf = as<arma::ivec>(parm["max_inf"]);
  parms.lamda_ex = as<arma::vec>(parm["lamda_ex"]);
  parms.times = as<arma::vec>(parm["times"]);
  
  parms.mg_actions = as<arma::uvec>(parm["mg_actions"]);
  parms.mg_actionlist = as<arma::uvec>(parm["mg_actionlist"]);
  parms.mg_levels = as<arma::vec>(parm["mg_levels"]);
  parms.mg_resprout = as<LogicalVector>(parm["mg_resprout"]);
  parms.mg_stages = as<arma::umat>(parm["mg_stages"]);
  
  int startcount = max(max(parms.K, parms.N), static_cast<int>(std::ceil(parms.K / arma::min(parms.omega(span(1, parms.omega.size() - 1))))));
  
  double time_max = parms.times(parms.times.n_elem - 1);
  //define functions based on options
  int dispersalfn = as<int>(parm["dispersalfn"]);
  int seedshadow = as<int>(parm["seedshadow"]);
  int beta_meth = as<int>(parm["beta_meth"]);
  
  if(dispersalfn == 0) {
    parms.m.fill(1/((parms.bbox(1) - parms.bbox(0))*(parms.bbox(3) - parms.bbox(2))));
    kernel1 = &flatdisp;
    kernel2 = &flatdisp2;
  } else if(dispersalfn == 1) {
    kernel1 = &expdisp;
    kernel2 = &expdisp2;
  } else if(dispersalfn == 2) {
    kernel1 = &fatdisp;
    kernel2 = &fatdisp2;
  } else if(dispersalfn == 3) {
    kernel1 = &normdisp;
    kernel2 = &normdisp2;
  }
  
  if(seedshadow == 0) {
    seedkern = &flatseed;
  } else if(seedshadow == 1) {
    seedkern = &expseed;;
  } else if(seedshadow == 3) {
    seedkern = &normseed;
  }
  
  if (beta_meth == 0) {
    beta_f = &beta_flat;
  } else if(beta_meth == 1) {
    beta_f = &beta_step;
  } else if(beta_meth == 2) {
    beta_f = &beta_lin;
  }

  //Bookkeeping values
	statelist state;
	state.treecount = parms.N;
  state.treeindex = parms.N - 1;
  state.IDcount = parms.N;
  state.ID = as<arma::uvec>(init["ID"]);
  state.X = as<arma::vec>(init["X"]);
  state.Y = as<arma::vec>(init["Y"]);
  state.S = as<arma::uvec>(init["Stage"]);
  state.I = as<arma::ivec>(init["Infections"]);
  state.F.zeros(startcount);
  state.B.zeros(startcount);
  state.time = parms.times(0);
  state.next_record = parms.times.begin();
  state.next_act = parms.mg_actions.begin();

  
  for(uword k = 0; k < state.treecount; k++) {
    state.B(k) = beta_f(parms.beta(state.S(k)), state.I(k), parms.max_inf(state.S(k)));
    state.F(k) = sum((kernel1(distance(state, k), parms.m(state.S)) % state.I % parms.lamda(state.S)));
  }

  state.E = fmax(0, 1 - sum(parms.omega(state.S)) / parms.K);
  arma::mat R(startcount, 6);
  recalc_probs(state, parms, R);
  
  RNGScope scope;
  IntegerVector Index = seq_len(startcount) - 1;
  IntegerVector actions = seq_len(6);
  uword s;
  int i;
  int j;
  int action;
  
  arma::mat printmatrix(startcount,6);  
  std::ofstream outfile;
  outfile.open(filename.c_str(), ios::out | ios::app);
  arma::rowvec jrow(6);
  std::ofstream outdiag;


//Start the loop
  while (*(state.next_record) < time_max) {

    //Calculate individual-level event rates and next time step.

    state.time = state.time + as<double>(rexp(1, accu(R)));

    //Record when we pass a value in the times vector
    while (state.time > *(state.next_record) && *(state.next_record) < time_max) {
      mgmt_action(state, parms, R);
      print_state(state, outfile, printmatrix);
      if (progress) {
       Rcpp::Rcout << "\nTime: " << *(state.next_record) << ", Population:" << state.treecount << ", Complete:" << Rf_fround(100 * *(state.next_record)/time_max, 1) << " %     ";
      }
     ++(state.next_record);
    }

    if(state.treecount == 0) break;
    //Select the individual that will change this time step
    j = as<uword>(Rcpp::RcppArmadillo::sample(Index, 1, false, as<NumericVector>(wrap(sum(R,1)))));
    jrow = R.row(j);
    action = as<int>(Rcpp::RcppArmadillo::sample(actions, 1, false, as<NumericVector>(wrap(jrow))));
    s = state.S(j);
    i = state.I(j);
    
    try {
    //Call the updating function on the individual
    switch (action) {  
    case 1: 
      reproduce(state, parms, R, j, s, state.treecount);
      break; 
    case 2:
      die(state, parms, R, j, s, state.treeindex);
      break;
    case 3: //grow
       grow(state, parms, R, j, s, i);
       break;
    case 4: //sicken
       sicken(state, parms, R, j, s, i);
       break;
    case 5: //recover
       recover(state, parms, R, j, s, i);
       break;
    case 6: //resprout
       resprout(state, parms, R, j, s, i);
       break;
    }
    
    } catch (std::exception &e) {
      if(diagnostics) {
      outdiag.open(diagfile.c_str(), ios::out | ios::app);

      R.save(outdiag, arma::csv_ascii);
      outdiag << "\n" << action << " " << j << " " << s << " " << i << " " << state.treeindex << " " << state.treecount << " " << state.F(j) <<
                  " " << state.B(j) << std::endl;
      outdiag << jrow << std::endl;
      outdiag << startcount << std::endl;
      outdiag << e.what() << std::endl;
      Rcerr << "Error in simulation loop. Diagnostics in " << diagfile << std::endl;
      outdiag.close();

      break;
      } else {
        Rcerr << "Error in simulation loop.";
        break;
      }
    };
    

  
  }
  if(state.treecount > 0) {
    print_state(state, outfile, printmatrix);
  }
  if (progress) {
     Rcpp::Rcout << "\nTime: " << *(state.next_record) << ", Population:" << state.treecount << ", Complete:" << Rf_fround(100 * *(state.next_record)/time_max, 1) << " %     ";
  }
  #if PROFILE
  ProfilerStop();
  #endif
  outfile.close();
  return 1;
}

void mgmt_action(statelist &state, parmlist &parms, arma::mat R) {
  uword act = parms.mg_actionlist(*(state.next_act));
  switch (act) {
    case 0:
      break;
    case 1:
      thin_evenly(state, parms, R, parms.mg_stages.row(act), parms.mg_levels(act), parms.mg_resprout(act));
      break;
    case 2:
      thin_spacing(state, parms, R, parms.mg_stages.row(act), parms.mg_levels(act), parms.mg_resprout(act));
      break;
  };
  
  ++(state.next_act);
}
