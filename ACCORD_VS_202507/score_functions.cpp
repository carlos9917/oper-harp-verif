#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double score_function_cpp(double phi, double ob, double k = 0.1, double A = 4.0) {
  if (NumericVector::is_na(phi) || NumericVector::is_na(ob)) {
    return NA_REAL;
  }

  if (ob > k) {
    if (phi < ob - k) {
      return phi / (ob - k);  // Equation (2a)
    } else if (phi <= ob) {
      return 1.0;  // Equation (2b)
    } else {  // phi > ob
      return std::max(1.0 - (phi - ob) / (A * ob), 0.0);  // Equation (2c)
    }
  } else {  // ob <= k
    if (phi <= k) {
      return 1.0;  // Equation (3a)
    } else {  // phi > k
      return std::max(1.0 - (phi - k) / (A * k), 0.0);  // Equation (3b)
    }
  }
}
