#define PROFILE 0
#define PROFILE_FILE "~/out.prof"

#if PROFILE
#include <gperftools/profiler.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


#include "data_structures.cpp"
#include "dispersal_functions"
#include "updating_functions"
#include "infection_density_dependence"

//' @export
// [[Rcpp::export]]
int run_sodi_rcpp(DataFrame init, List parms, bool progress) {
  
  #if PROFILE
  ProfilerStart(PROFILE_FILE);
  #endif
  
  //populate the parameter structure from the input list
  parmlist parms;
  parms.K = as<int>(parms["K"]);
  parms.N = as<int>(parms["n0"]);
  parm.bbox = as<NumericVector>(parms["bbox"]);
  parms.seedm = as<arma::vec>(parms["seedm"]);
  parms.m = as<arma::vec>(parms["m"]);
  parms.f = as<arma::vec>(parms["f"]);
  parms.g = as<arma::vec>(parms["g"]);
  parms.d = as<arma::vec>(parms["d"]);
  parms.r = as<arma::vec>(parms["r"]);
  parms.alpha = as<arma::vec>(parms["alpha"]);
  parms.lamda = as<arma::vec>(parms["lamda"]);
  parms.beta = as<arma::vec>(parms["beta"]);
  parms.mu = as<arma::vec>(parms["mu"]);
  parms.xi = as<arma::vec>(parms["xi"]);
  parms.omega = as<arma::vec>(parms["omega"]);
  parms.ss = as<arma::uvec>(parms["ss"]);
  parms.max_inf = as<arma::ivec>(parms["max_inf"]);


  //define functions based on options
  int dispersalfn = as<int>(parms["dispersalfn"]);
  int seedshadow = as<int>(parms["seedshadow"]);
  int beta_meth = as<int>(parms["beta_meth"]);
  
  arma::vec (*kernel)(arma::vec d, arma::vec m);
  arma::vec (*kernel2)(arma::vec d, double m);
  if(dispersalfn == 0) {
    m.fill(1/((bbox(1) - bbox(0))*(bbox(3) - bbox(2))));
    kernel = flatdisp;
    kernel2 = flatdisp2;
  } else if(dispersalfn == 1) {
    kernel = expdisp;
    kernel2 = expdisp2;
  } else if(dispersalfn == 2) {
    kernel = fatdisp;
    kernel2 = fatdisp2;
  } else if(dispersalfn == 3) {
    kernel = normdisp;
    kernel2 = normdisp2;
  }
  
  arma::vec (*seedkern)(double x, double y, int s, parmlist parms);
  if(seedshadow == 0) {
    seedkern = flatseed;
  if(seedshadow == 1) {
    seedkern = expseed;;
  } else if(seedshadow == 3) {
    seedkern = normseed;
  }
  
  double (*beta_f)(double beta, int i, int max_inf);
  if (beta_meth == 0) {
    beta_f = beta_flat;
  } else if(beta_meth == 1) {
    beta_f = beta_step;
  } else if(beta_meth == 2) {
    beta_f = beta_lin;
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
  state.F.zeros(K);
  state.B.zeros(K);
  state.time = times(0);
  
  for(uword k = 0; k < state.treecount; k++) {
    state.B(k) = beta_f(parms.beta(state.S(k)), state.I(k), parms.max_inf(state.S(k)));
    state.F(k) = sum(kernel(distance(state.X, state.Y, k), parms.m(state.S)) % state.I % parms.lamda(state.S)) * state.B(k);
  }
  
  state.E = fmax(0, 1 - sum(parms.omega(state.S)) / parms.K);
  arma::mat R(parms.K, 6);
  R.col(1) = state.E * parms.f(state.S);
  R.col(2) = parms.d(state.S) + state.I % alpha(state.S) % (1 - parms.r(state.S));
  R.col(3) = parms.g(state.S);
  R.col(4) = state.F;
  R.col(5) = state.I % parms.mu(state.S);
  R.col(6) = state.I % alpha(state.S) % parms.r(state.S));
  
  
  NumericVector times = as<NumericVector>(parms["times"]);
  CharacterVector timenames = as<CharacterVector>(parms["timenames"]);
  
  double time_max = times(times.length() - 1);
  NumericVector::iterator next_record = times.begin();
  CharacterVector::iterator next_timename = timenames.begin();

  
  //Interim values
  RNGScope scope;
  uword s;
  int i;
  int j;
  int action;
  double theta;
  
  
  //Assign columns of the data frame to vectors of state value




  //Initiate datafrae

  //Initiate iterators
  ++next_record;
  ++next_timename;

  mat printmatrix(K,6);  
  ofstream outfile;
  outfile.open(filename, ios::out | ios::app);
//Start the loop
  while (time < time_max) {
 
    //Record when we pass a value in the times vector
    if (time > *next_record) {
      print_state(*next_record, state, outfile, printmatrix)
      void print_state(double outtime, statelist &state, ofstream &outfile, mat &printmatrix) {
        printmartix.col(0) = outtime;
        printmatrix.col(1) = state.ID;
        printmatrix.col(2) = state.X;
        printmatrix.col(3) = state.Y;
        printmatrix.col(4) = state.S;
        printmatrix.col(5) = state.I;
        printmatrix(span(0, state.treeindex), span(0, 5)).save(outfile, csv_ascii);
      }  
      
      ++next_record;
      ++next_timename;

      }

    //Calculate individual-level event rates and next time step.
    
    time = time + as<double>(rexp(1, accu(R)));
    //T.fill(time);
    
    
    //Select the individual that will change this time step
    j = as<uword>(Rcpp::RcppArmadillo::sample(Index, 1, false, as<NumericVector>(wrap(sum(state.R, dim=1)))));
    action = as<int>(Rcpp::RcppArmadillo::sample(actions, 1, false, as<NumericVector(wrap(state.R.row(j)))));
    s = state.S(j);
    i = state.I(j);

    switch (action) {  
    case 1: 
      reproduce(&state, &parms, &R, j, s, state.treecount);
      break; 
    case 2:
    die(&state, &parms, &R, j, s, state.treeindex);
    break;
    case 3: //grow
       grow(&state, &parms, &R, j, s, i);
       break;
    case 4: //sicken
       sicken(&state, &parms, &R, j, s, i);
       break;

    case 5: //recover
       recover(&state, &parms, &R, j, s, i);
       break;
    case 6: //resprout
       resprout(&state, &parms, &R, j, s, i);
       break;
    }

  F = as<arma::vec>(wrap(pmax(0, as<NumericVector>(wrap(F)))));
  
  }
  
  #if PROFILE
  ProfilerStop();
  #endif
  outfile.close();
  return 1;
}

