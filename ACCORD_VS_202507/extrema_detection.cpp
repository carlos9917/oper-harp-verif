#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix find_local_extrema_cpp(NumericMatrix arr, String mode, double tolerance = 0.0) {
  int rows = arr.nrow();
  int cols = arr.ncol();

  // Use vector to store results, then convert to matrix
  std::vector<double> extrema_rows, extrema_cols, extrema_vals;

  // Check each point against its 3x3 neighborhood (excluding boundaries)
  for (int i = 1; i < rows - 1; i++) {
    for (int j = 1; j < cols - 1; j++) {
      double current_val = arr(i, j);

      // Skip NA values
      if (NumericVector::is_na(current_val)) continue;

      bool is_extremum = false;

      if (mode == "max") {
        // Local maximum: current value >= all neighbors within tolerance
        double max_neighbor = current_val;
        for (int di = -1; di <= 1; di++) {
          for (int dj = -1; dj <= 1; dj++) {
            double neighbor_val = arr(i + di, j + dj);
            if (!NumericVector::is_na(neighbor_val) && neighbor_val > max_neighbor) {
              max_neighbor = neighbor_val;
            }
          }
        }
        is_extremum = (current_val >= max_neighbor - tolerance) && (current_val == max_neighbor);
      } else {
        // Local minimum: current value <= all neighbors within tolerance
        double min_neighbor = current_val;
        for (int di = -1; di <= 1; di++) {
          for (int dj = -1; dj <= 1; dj++) {
            double neighbor_val = arr(i + di, j + dj);
            if (!NumericVector::is_na(neighbor_val) && neighbor_val < min_neighbor) {
              min_neighbor = neighbor_val;
            }
          }
        }
        is_extremum = (current_val <= min_neighbor + tolerance) && (current_val == min_neighbor);
      }

      if (is_extremum) {
        extrema_rows.push_back(i + 1);  // Convert to R indexing (1-based)
        extrema_cols.push_back(j + 1);
        extrema_vals.push_back(current_val);
      }
    }
  }

  // Convert to matrix
  int n_extrema = extrema_rows.size();
  NumericMatrix result(n_extrema, 3);
  colnames(result) = CharacterVector::create("row", "col", "value");

  for (int k = 0; k < n_extrema; k++) {
    result(k, 0) = extrema_rows[k];
    result(k, 1) = extrema_cols[k];
    result(k, 2) = extrema_vals[k];
  }

  return result;
}
