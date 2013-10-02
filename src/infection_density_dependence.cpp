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
  
double beta_lin(double beta, int i, int max_inf) {
  if(i >= max_inf) {
    return 0;
  } else {
    return beta * (1 - i / max_inf);
  }
}
