#include "includes.h"
#include "data_structures.h"
using namespace Rcpp;


double lamda_interp(statelist &state, parmlist &parms) {
  
  arma::uword index = state.next_record - parms.times.begin();
  double a = parms.lamda_ex(index-1) + (state.time - parms.times(index-1)) * (parms.lamda_ex(index) - parms.lamda_ex(index-1)) / (parms.times(index) - parms.times(index-1));
//  Rcout << "\n" << index << " " << parms.times(index-1) << " " << parms.times(index) << " " << state.time << " " << parms.lamda_ex(index-1) << " " << parms.lamda_ex(index) << " " << a;
  return a;
}

//Lamda-ex interpolate Delta
//Function of state.time and parms (specifally times and lamda_x)
//lamda ex must either have length one or length = length(times)
//Check at beginning of calculation, define function depending on length
//Save last time point, calculate delta, add it to the Force of infection column
//Calculate some equivalent tree-distances


//Steps
// Determine current time
// Determine last recorded time and next recorded time
// Determine value at last and next recorded time
// value = val(before) + val(time)*((val(after) - val(before))/(time(after) - time(before)))
