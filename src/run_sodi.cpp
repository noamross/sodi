// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

//The dispersal functions 

//' @export
// [[Rcpp::export]]
arma::vec flatdisp(arma::vec d, arma::vec m) {
  return m;
}

//' @export
// [[Rcpp::export]]
arma::vec flatdisp2(arma::vec d, double m) {
  d.fill(m);
  return d;
}

//' @export
// [[Rcpp::export]]
arma::vec expdisp(arma::vec d, arma::vec m) {
  return (square(m) / (2*M_PI)) % exp(-m % d);
}

//' @export
// [[Rcpp::export]]
arma::vec expdisp2(arma::vec d, double m) {
  return (m*m / (2*M_PI)) * exp(-m * d);
}

//' @export
// [[Rcpp::export]]
arma::vec normdisp(arma::vec d, arma::vec m) {
  return (square(m) / sqrt(2*M_PI)) % exp(-0.5 * square(m % d));
}

//' @export
// [[Rcpp::export]]
arma::vec normdisp2(arma::vec d, double m) {
  return (m*m / sqrt(2*M_PI)) * exp(-0.5 * square(m * d));
}

//' @export
// [[Rcpp::export]]
arma::vec fatdisp(arma::vec d, arma::vec m) {
  return (square(m) / (24* M_PI)) % exp(- sqrt(m % d));
}

//' @export
// [[Rcpp::export]]
arma::vec fatdisp2(arma::vec d, double m) {
  return (m*m / (24* M_PI)) * exp(- sqrt(m * d));
}

arma::vec distance(arma::vec X, arma::vec Y, int j) {
  return sqrt(square(X - X(j)) + square(Y - Y(j)));
}

arma::vec powvec(arma::vec A, arma::ivec p) {
  uword n = A.n_elem;
  arma::vec B(n);
  for(uword i = 0; i < n; ++i) {
    B(i) = pow(A(i), p(i));
  }
  return B;
}

double rnorm1(double p) {
  return as<double>(rnorm(1, 0, p));
}

double rexp1(double p) {
  return as<double>(rexp(1, p));
}

//' @export
// [[Rcpp::export]]
List run_sodi_rcpp(DataFrame init, List parms, bool progress) {
  //Initiate the outputlist
  //Unpack the parameters

  int K = as<int>(parms["K"]);
  int N = as<int>(parms["n0"]);
  int max_inf = as<int>(parms["max_inf"]);
  uword treeindex = N - 1;
  uword treecount = N;
  uword IDcount = N;
  NumericVector bbox = as<NumericVector>(parms["bbox"]);
  int dispersalfn = as<int>(parms["dispersalfn"]);
  int seedshadow = as<int>(parms["dispersalfn"]);
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
  double theta;
  
  RNGScope scope;
  double E;
  uword s;
  int i;
  int j;
  NumericVector times = as<NumericVector>(parms["times"]);
  CharacterVector timenames = as<CharacterVector>(parms["timenames"]);
  double time = times(0);
  double time_max = times(times.length() - 1);
  NumericVector::iterator next_record = times.begin();
  CharacterVector::iterator next_timename = timenames.begin();

  
  arma::vec irates;
  
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
  
  //Assign columns of the data frame a vectors.
  arma::uvec ID = as<arma::uvec>(init["ID"]);
  arma::vec X = as<arma::vec>(init["X"]);
  arma::vec Y = as<arma::vec>(init["Y"]);
  arma::uvec S = as<arma::uvec>(init["Stage"]);
  arma::ivec I = as<arma::ivec>(init["Infections"]);
  arma::vec F(K);
  F.fill(0);
  for(uword k = 0; k < treecount; k++) {
      F(k) = sum(kernel(distance(X, Y, k), m(S)) % I % lamda(S)) * beta(S(k));
  }
  arma::vec T(K);
  T.fill(time);
  init = Rcpp::DataFrame::create(Named("Time", T), init, Named("Force", F));
  List sodi = List::create(Named("0", init));

  IntegerVector Index = seq_len(K) - 1;
  IntegerVector actions = seq_len(6);
  //unpack the parms list so all the variables are scalars/vectors
  //while the time is less than the last value of the times sequence
  
  
 // Progress p(times.length(), progress);
  ++next_record;
  ++next_timename;
  
  //Initiate a progress bar
//  if (progress) {
//  Environment e = Environment::global_env(); 
//  Function sodistep = e["sodistep"];
//  }
  
  
  while (time < time_max) {
//    if (Progress::check_abort())
  //    return;
    //Record when we pass a value in the times vector
    
    if (time > *next_record) {
    sodi.push_back(Rcpp::DataFrame::create(_["Time"]=T,
                                           _["ID"]=ID,
                                           _["X"]=X,
                                           _["Y"]=Y,
                                           _["Stage"]=S,
                                           _["Infections"]=I,
                                           _["Force"]=F));
    ++next_record;
    ++next_timename;
    if (progress) {
      Rf_PrintValue(wrap(time/time_max));
//      sodistep();
    }

//    p.increment();
    }

    
    double E = 1 - sum(omega(S)) / K;
    irates = E * f(S) % powvec(xi(S), I) + d(S) + g(S) + F + I % mu(S) + I % alpha(S);
    time = time + as<double>(rexp(1, sum(irates)));
    T.fill(time);

    uword j = as<uword>(Rcpp::RcppArmadillo::sample(Index, 1, false, as<NumericVector>(wrap(irates))));
    s = S(j);
    i = I(j);




    NumericVector possibs = NumericVector::create(
      E * f(s) * pow(xi(s), i), //1. reproduce
      d(s) + i * alpha(s) * (1 - r(s)), //2. die/succumb
      g(s), //3. grow
      F(j), //4. sicken
      i * mu(s), //5. recover
      i * alpha(s) * r(s)  //6. succumb + resprout
    );



    
    switch (as<int>(Rcpp::RcppArmadillo::sample(actions, 1, false, possibs))) {
    case 1: //reproduce
       IDcount += 1;
       if (seedshadow == 0) {
        X(treecount) = as<double>(runif(1, bbox(0), bbox(1)));
        Y(treecount) = as<double>(runif(1, bbox(2), bbox(3)));
       } else if(seedshadow > 0) {
         theta = as<double>(runif(1, 0, 2*M_PI));
         seeddist = seedkern(seedm(s));
         X(treecount) = X(j) + cos(theta) * seeddist;
         Y(treecount) = Y(j) + sin(theta) * seeddist;
       }
       S(treecount) = 1;
       F(treecount) = sum(kernel(distance(X, Y, treecount), m(S)) % I % lamda(S)) * beta(s) ;
       ID(treecount) = IDcount;
       treecount += 1;
       treeindex += 1;
       break; 
    case 2: //die
       F -= kernel2(distance(X, Y, j), m(s)) % beta(S) * i * lamda(s);
       ID(j) = ID(treeindex);
       X(j) = X(treeindex);
       Y(j) = Y(treeindex);
       S(j) = S(treeindex);
       I(j) = I(treeindex);
       F(j) = F(treeindex);
       ID(treeindex) = 0;
       X(treeindex) = 0;
       Y(treeindex) = 0;
       S(treeindex) = 0;
       I(treeindex) = 0;
       F(treeindex) = 0;
       treecount -= 1;
       treeindex -= 1;
       break;
    case 3: //grow
       F += kernel2(distance(X, Y, j), m(s)) % beta(S) * i * (lamda(s+1) - lamda(s));
       S(j) += 1;
       break;
    case 4: //sicken
       if (I(j) < max_inf) {
         I(j) += 1;
         F += kernel2(distance(X, Y, j), m(s)) % beta(S) * lamda(s);
       }
       break;
    case 5: //recover
       I(j) -= 1;
       F -= kernel2(distance(X, Y, j), m(s)) % beta(S) * lamda(s);
       break;
    case 6: //resprout
       F -= kernel2(distance(X, Y, j), m(s)) % beta(S) * i * lamda(s);
       I(j) = 0;
       S(j) = ss(s);  //Set the stage back to the "base stage" for the species
       
    }
    
    //RcppProgressBar
  }
  
  return sodi;
}

//post processing
//rbind the list of matrices, eliminate zero rows
//Create records for dead trees?
//Convert state to stage + species factors


