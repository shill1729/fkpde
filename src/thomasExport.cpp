#include <Rcpp.h>
#include "Tridiag.h"
using namespace Rcpp;


//' Thomas algorithm for solving tridiagonal linear systems. Stability only for diagonally dominant or symmetric positive definite.
//'
//' @param a lower diagonal of matrix
//' @param b diagonal of matrix
//' @param c upper diagonal of matrix
//' @param d vector on RHS of equation
//'
//' @description {Thomas algorithm for tridiagonal linear systems. This can be used to solve
//' PDEs arising from Feynman-Kac connections.}
//' @details {The algorithm is on wikipedia and quite straightforward.}
// [[Rcpp::export]]
std::vector<double> thomasAlgorithm(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d)
{
  return Tridiag::thomasAlgorithm(a, b, c, d);
}
