#include <iostream>
#include <fstream>
#include <RcppArmadillo.h>
#include "data_structures.h"
#include "print_state.h"

void print_state(double outtime, statelist &state, std::ofstream outfile, arma::mat &printmatrix) {
  printmatrix.col(0).fill(outtime);
  printmatrix.col(1) = arma::conv_to<arma::vec>::from(state.ID);
  printmatrix.col(2) = state.X;
  printmatrix.col(3) = state.Y;
  printmatrix.col(4) = arma::conv_to<arma::vec>::from(state.S);
  printmatrix.col(5) = arma::conv_to<arma::vec>::from(state.I);
  arma::mat p = printmatrix(arma::span(0, state.treeindex), arma::span(0, 5));
  p.save(outfile, arma::csv_ascii);
}  