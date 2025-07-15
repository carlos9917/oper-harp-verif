#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double get_neighbourhood_extreme_cpp(NumericMatrix arr, int i, int j, int L, String mode) {
  int rows = arr.nrow();
  int cols = arr.ncol();

  // Convert to 0-based indexing
  i = i - 1;
  j = j - 1;

  // Define neighborhood bounds
  int i_min = std::max(0, i - L);
  int i_max = std::min(rows - 1, i + L);
  int j_min = std::max(0, j - L);
  int j_max = std::min(cols - 1, j + L);

  double extreme_val;
  bool first_valid = true;

  for (int row = i_min; row <= i_max; row++) {
    for (int col = j_min; col <= j_max; col++) {
      double val = arr(row, col);
      if (!NumericVector::is_na(val)) {
        if (first_valid) {
          extreme_val = val;
          first_valid = false;
        } else {
          if (mode == "max" && val > extreme_val) {
            extreme_val = val;
          } else if (mode == "min" && val < extreme_val) {
            extreme_val = val;
          }
        }
      }
    }
  }

  return extreme_val;
}
