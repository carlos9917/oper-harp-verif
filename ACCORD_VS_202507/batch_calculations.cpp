#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calculate_component_scores_cpp(NumericMatrix extrema, NumericMatrix field, 
                                           int L, String mode, String component_type,
                                           double k = 0.1, double A = 4.0) {
  int n_extrema = extrema.nrow();
  NumericVector scores(n_extrema);

  for (int idx = 0; idx < n_extrema; idx++) {
    int i = extrema(idx, 0);  // Already 1-based from R
    int j = extrema(idx, 1);
    double extrema_val = extrema(idx, 2);

    // Get neighborhood extreme using the same logic as get_neighbourhood_extreme_cpp
    int rows = field.nrow();
    int cols = field.ncol();

    // Convert to 0-based indexing
    int i_zero = i - 1;
    int j_zero = j - 1;

    // Define neighborhood bounds
    int i_min = std::max(0, i_zero - L);
    int i_max = std::min(rows - 1, i_zero + L);
    int j_min = std::max(0, j_zero - L);
    int j_max = std::min(cols - 1, j_zero + L);

    double field_extreme;
    bool first_valid = true;

    for (int row = i_min; row <= i_max; row++) {
      for (int col = j_min; col <= j_max; col++) {
        double val = field(row, col);
        if (!NumericVector::is_na(val)) {
          if (first_valid) {
            field_extreme = val;
            first_valid = false;
          } else {
            if (mode == "max" && val > field_extreme) {
              field_extreme = val;
            } else if (mode == "min" && val < field_extreme) {
              field_extreme = val;
            }
          }
        }
      }
    }

    // Calculate score using the same logic as score_function_cpp
    double score;
    double phi, ob;

    if (component_type == "obs_based") {
      // For SLX_ob_max and SLX_ob_min: score forecast against observed extrema
      phi = field_extreme;
      ob = extrema_val;
    } else {
      // For SLX_fc_max and SLX_fc_min: score forecast extrema against observed field
      phi = extrema_val;
      ob = field_extreme;
    }

    if (NumericVector::is_na(phi) || NumericVector::is_na(ob)) {
      score = NA_REAL;
    } else if (ob > k) {
      if (phi < ob - k) {
        score = phi / (ob - k);  // Equation (2a)
      } else if (phi <= ob) {
        score = 1.0;  // Equation (2b)
      } else {  // phi > ob
        score = std::max(1.0 - (phi - ob) / (A * ob), 0.0);  // Equation (2c)
      }
    } else {  // ob <= k
      if (phi <= k) {
        score = 1.0;  // Equation (3a)
      } else {  // phi > k
        score = std::max(1.0 - (phi - k) / (A * k), 0.0);  // Equation (3b)
      }
    }

    scores[idx] = score;
  }

  return scores;
}
