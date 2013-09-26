// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
//#include <gperftools/profiler.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

//The dispersal functions. 

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


//Distance calculation.  Returns vector of Euclidian distances from element j
arma::vec distance(arma::vec X, arma::vec Y, int j) {
  return sqrt(square(X - X(j)) + square(Y - Y(j)));
}

//Vectorized power calculation
arma::vec powvec(arma::vec A, arma::ivec p) {
  uword n = A.n_elem;
  arma::vec B(n);
  for(uword i = 0; i < n; ++i) {
    B(i) = pow(A(i), p(i));
  }
  return B;
}

//Random number generators
double rnorm1(double p) {
  return as<double>(rnorm(1, 0, p));
}

double rexp1(double p) {
  return as<double>(rexp(1, p));
}


//Functions to determine infection probability based on infection number
double beta_flat(double beta, int i, int max_inf) {
  return beta;
}

double beta_step(double beta, int i, int max_inf) {
  if(i >= max_inf) {
    return 0;
  } else {
    return beta;
  }
}
  
double beta_lin(double beta, int i, int max_inf)  {
  if(i >= max_inf) {
    return 0;
  } else {
    return beta * (1 - i / max_inf);
  }
}


//' @export
// [[Rcpp::export]]
List run_sodi_rcpp(DataFrame init, List parms, bool progress) {
  //ProfilerStart("~/code/sodi/run_sodi.prof");
  //Unpack the parameters

  int K = as<int>(parms["K"]);
  int N = as<int>(parms["n0"]);
  NumericVector bbox = as<NumericVector>(parms["bbox"]);

  //stage-specific parameters
  arma::vec seedm = as<arma::vec>(parms["seedm"]);
  arma::vec m = as<arma::vec>(parms["m"]);
  arma::vec f = as<arma::vec>(parms["f"]);
  arma::vec g = as<arma::vec>(parms["g"]);
  arma::vec d = as<arma::vec>(parms["d"]);
  arma::vec r = as<arma::vec>(parms["r"]);
  arma::vec alpha = as<arma::vec>(parms["alpha"]);
  arma::vec lamda = as<arma::vec>(parms["lamda"]);
  arma::vec beta = as<arma::vec>(parms["beta"]);
  arma::vec mu = as<arma::vec>(parms["mu"]);
  arma::vec xi = as<arma::vec>(parms["xi"]);
  arma::vec omega = as<arma::vec>(parms["omega"]);
  arma::uvec ss = as<arma::uvec>(parms["ss"]);
  arma::ivec max_inf = as<arma::ivec>(parms["max_inf"]);

  //Options to select functions for dispersal and infection
  int dispersalfn = as<int>(parms["dispersalfn"]);
  int seedshadow = as<int>(parms["seedshadow"]);
  int beta_meth = as<int>(parms["beta_meth"]);
  
  //define functions based on options
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
  
  double (*seedkern)(double p);
  if(seedshadow == 1) {
    seedkern = rexp1;
  } else if(seedshadow == 3) {
    seedkern = rnorm1;
  }
  double seeddist;
  
  double (*beta_f)(double beta, int i, int max_inf);
  if (beta_meth == 0) {
    beta_f = beta_flat;
  } else if(beta_meth == 1) {
    beta_f = beta_step;
  } else if(beta_meth == 2) {
    beta_f = beta_lin;
  }

  //Bookkeeping values
  uword treeindex = N - 1;
  uword treecount = N;
  uword IDcount = N;
  NumericVector times = as<NumericVector>(parms["times"]);
  CharacterVector timenames = as<CharacterVector>(parms["timenames"]);
  double time = times(0);
  double time_max = times(times.length() - 1);
  NumericVector::iterator next_record = times.begin();
  CharacterVector::iterator next_timename = timenames.begin();
  arma::vec T(K);
  T.fill(time);
  IntegerVector Index = seq_len(K) - 1;
  IntegerVector actions = seq_len(6);
  int action;
  arma::rowvec status;
  
  //Interim values
  RNGScope scope;
  double E;
  uword s;
  int i;
  int j;
  double theta;
  
  
  //Assign columns of the data frame to vectors of state value
  arma::uvec ID = as<arma::uvec>(init["ID"]);
  arma::vec X = as<arma::vec>(init["X"]);
  arma::vec Y = as<arma::vec>(init["Y"]);
  arma::uvec S = as<arma::uvec>(init["Stage"]);
  arma::ivec I = as<arma::ivec>(init["Infections"]);
  
  //Intermediate state values 
  arma::vec irates(K);
  arma::vec F(K);
  arma::vec beta_i(K);
  irates.fill(0);
  F.fill(0);
  beta_i.fill(0);

  for(uword k = 0; k < treecount; k++) {
      beta_i(k) = beta_f(beta(S(k)), I(k), max_inf(S(k)));
      F(k) = sum(kernel(distance(X, Y, k), m(S)) % I % lamda(S)) * beta_i(k);
  }
  

  //Initiate datafrae
  init = Rcpp::DataFrame::create(Named("Time", T), init,  Named("Force", F));
  List sodi = List::create(Named("0", init));


  //Initiate iterators
  ++next_record;
  ++next_timename;
  
//Start the loop
  while (time < time_max) {
 
    //Record when we pass a value in the times vector
    if (time > *next_record) {
      sodi.push_back(Rcpp::DataFrame::create(_["Time"]=wrap(T.fill(*next_record)),
                                             _["ID"]=wrap(ID),
                                             _["X"]=wrap(X),
                                             _["Y"]=wrap(Y),
                                             _["Stage"]=wrap(S),
                                             _["Infections"]=wrap(I),
                                             //_["Beta"]=beta_i,
                                             _["Force"]=wrap(F)));
      ++next_record;
      ++next_timename;
      if (progress) {
        Rf_PrintValue(wrap(time/time_max));
      }
    }

    //Calculate individual-level event rates and next time step.
    double E = fmax(0, 1 - sum(omega(S)) / K);
    irates = E * f(S) % powvec(xi(S), I) + d(S) + g(S) + F + I % mu(S) + I % alpha(S);
    time = time + as<double>(rexp(1, sum(irates)));
    //T.fill(time);
    
    
    //Select the individual that will change this time step
    uword j = as<uword>(Rcpp::RcppArmadillo::sample(Index, 1, false, as<NumericVector>(wrap(irates))));
    s = S(j);
    i = I(j);

    //Calcualate the probabilities of different events for the individual
    NumericVector possibs = NumericVector::create(
      E * f(s) * pow(xi(s), i), //1. reproduce
      d(s) + i * alpha(s) * (1 - r(s)), //2. die/succumb
      g(s), //3. grow
      F(j), //4. sicken
      i * mu(s), //5. recover
      i * alpha(s) * r(s)  //6. succumb + resprout
    );
    
    
    //Determine the event
    action = as<int>(Rcpp::RcppArmadillo::sample(actions, 1, false, possibs));
   // status << j << action << S(j) << I(j) << beta_i(j) << F(j);
   // Rf_PrintValue(wrap(status));

    
    switch (action) {
    case 1: //reproduce, create a new tree
      IDcount += 1;
      ID(treecount) = IDcount;
      if (seedshadow == 0) {
        X(treecount) = as<double>(runif(1, bbox(0), bbox(1)));
        Y(treecount) = as<double>(runif(1, bbox(2), bbox(3)));
      } else if(seedshadow > 0) {
        theta = as<double>(runif(1, 0, 2*M_PI));
        seeddist = seedkern(seedm(s));
        X(treecount) = X(j) + cos(theta) * seeddist;
        Y(treecount) = Y(j) + sin(theta) * seeddist;
      }
      S(treecount) = ss(s);
      I(treecount) = 0;
      beta_i(treecount) = beta(s);
      F(treecount) = sum(kernel(distance(X, Y, treecount), m(S)) % I % lamda(S)) * beta_i(treecount);
      treecount += 1;
      treeindex += 1;
      break; 
    case 2: //die
       F -= kernel2(distance(X, Y, j), m(s)) % beta_i * i * lamda(s);
       ID(j) = ID(treeindex);
       X(j) = X(treeindex);
       Y(j) = Y(treeindex);
       S(j) = S(treeindex);
       I(j) = I(treeindex);
       beta_i(j) = beta_i(treeindex);
       F(j) = F(treeindex);
       ID(treeindex) = 0;
       X(treeindex) = 0;
       Y(treeindex) = 0;
       S(treeindex) = 0;
       I(treeindex) = 0;
       beta_i(j) = 0;
       F(treeindex) = 0;
       treecount -= 1;
       treeindex -= 1;
       break;
    case 3: //grow
       S(j) += 1;
       beta_i(j) = beta_f(beta(S(j)), i, max_inf(S(j)));
       F(j) = F(j) * beta(S(j))/beta(s);
       F += kernel2(distance(X, Y, j), m(S(j))) % beta_i * i * (lamda(S(j)) - lamda(s));
       break;
    case 4: //sicken
       I(j) += 1;
       F += kernel2(distance(X, Y, j), m(s)) % beta_i * lamda(s);
       beta_i(j) = beta_f(beta(s), I(j), max_inf(s));
       F(j) = F(j) * beta_i(j) / beta_f(beta(s), i, max_inf(s));
       break;
    case 5: //recover
       I(j) -= 1;
       F -= kernel2(distance(X, Y, j), m(s)) % beta_i * lamda(s);
       beta_i(j) = beta_f(beta(s), I(j), max_inf(s));
       F(j) = sum(kernel(distance(X, Y, j), m(S)) % I % lamda(S)) * beta_i(j);
       break;
    case 6: //resprout
       S(j) = ss(s);  //Set the stage back to the "base stage" for the species
       I(j) = 0;
       beta_i(j) = beta((S(j)));
       F -= kernel2(distance(X, Y, j), m(s)) % beta_i * i * lamda(s);
       F(j) = sum(kernel(distance(X, Y, j), m(S)) % I % lamda(S)) * beta_i(j) ;
    }

  F = as<arma::vec>(wrap(pmax(0, as<NumericVector>(wrap(F)))));
  
  }
  //ProfilerStop();
  return sodi;

}

