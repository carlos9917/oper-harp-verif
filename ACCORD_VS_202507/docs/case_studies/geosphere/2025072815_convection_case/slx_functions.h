#ifndef SLX_FUNCTIONS_H
#define SLX_FUNCTIONS_H

#include <Rcpp.h>
using namespace Rcpp;

// Function declarations
double score_function_cpp_helper(double phi, double ob, double k = 0.1, double A = 4.0);
double get_neighbourhood_extreme_cpp_helper(NumericMatrix arr, int i, int j, int L, String mode);

#endif
