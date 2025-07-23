# Validation case for the SLX implementation
# Simple R implementation following the Python approach
# Based on test_simple_domain.py
# This one gets rid of the issues in the local search
# and works with the simple test case in Bent's paper


# SLX score implementation
S_score <- function(ob, phi, k = 0.1, A = 4) {
  # Piece-wise linear score function from Sass (2021)
  if (ob > k) {
    if (phi < ob - k) {
      return(phi / (ob - k))
    } else if (phi <= ob) {
      return(1.0)
    } else {
      return(max(1 - (phi - ob) / (A * ob), 0.0))
    }
  } else {  # ob <= k
    if (phi <= k) {
      return(1.0)
    } else {
      return(max(1 - (phi - k) / (A * k), 0.0))
    }
  }
}

#local_extreme_indices <- function(field, mode = 'max', delta = 0.0) {
#  # Return indices where the value is a global extreme (simplified approach), following python script
#  if (mode == 'max') {
#    target <- max(field, na.rm = TRUE)
#    mask <- abs(field - target) <= delta
#  } else {
#    target <- min(field, na.rm = TRUE)
#    mask <- abs(field - target) <= delta
#  }
#
#  # Get indices where mask is TRUE
#  indices <- which(mask, arr.ind = TRUE)
#  return(indices)
#}


local_extreme_indices <- function(field, mode = "max", tolerance = 0.0) {
#truly local version
# The one I was using before considers a point a local maximum if it is greater than or equal to all its neighbors.
# On a plateau (e.g., a large region where all values are 5), every point on the plateau will satisfy this condition, so all are marked as maxima.
# The version below considers the case to mark only the edge of a plateau as local maxima (i.e., only those points that are equal to their neighbors, but at least one neighbor is strictly less), you need to add an extra condition:
# The point is greater than or equal to all neighbors (as before)
# AND at least one neighbor is strictly less than the center value
# The original implementation, only marked the edge of the plateau (the "rim" of the block) as maxima, not the entire flat region.
  nrow_f <- nrow(field)
  ncol_f <- ncol(field)
  extrema <- matrix(numeric(0), 0, 3, dimnames = list(NULL, c("row", "col", "value")))

  for (i in 2:(nrow_f - 1)) {
    for (j in 2:(ncol_f - 1)) {
      val <- field[i, j]
      if (is.na(val)) next
      neighborhood <- field[(i-1):(i+1), (j-1):(j+1)]
      neighbors <- as.vector(neighborhood)
      neighbors <- neighbors[-5]  # Remove the center point

      if (mode == "max") {
        if (all(val >= neighbors - tolerance) && any(val > neighbors + tolerance)) {
          extrema <- rbind(extrema, c(i, j, val))
        }
      } else if (mode == "min") {
        if (all(val <= neighbors + tolerance) && any(val < neighbors - tolerance)) {
          extrema <- rbind(extrema, c(i, j, val))
        }
      }
    }
  }
  colnames(extrema) <- c("row", "col", "value")
  return(extrema)
}


#local_extreme_indices <- function(field, mode = "max", tolerance = 0.0) {
#  nrow_f <- nrow(field)
#  ncol_f <- ncol(field)
#  extrema <- matrix(numeric(0), 0, 3, dimnames = list(NULL, c("row", "col", "value")))
#
#  # Loop over interior points (avoid boundaries)
#  for (i in 2:(nrow_f - 1)) {
#    for (j in 2:(ncol_f - 1)) {
#      val <- field[i, j]
#      if (is.na(val)) next
#      neighborhood <- field[(i-1):(i+1), (j-1):(j+1)]
#      neighbors <- as.vector(neighborhood)
#      neighbors <- neighbors[-5]  # Remove the center point
#
#      if (mode == "max") {
#        if (all(val >= neighbors - tolerance)) {
#          extrema <- rbind(extrema, c(i, j, val))
#        }
#      } else if (mode == "min") {
#        if (all(val <= neighbors + tolerance)) {
#          extrema <- rbind(extrema, c(i, j, val))
#        }
#      }
#    }
#  }
#  colnames(extrema) <- c("row", "col", "value")
#  #colnames(extrema) <- c("row", "col")
#  return(extrema)
#}


neighbourhood_view <- function(field, i, j, L, mode = 'max') {
  # Get max/min value in (2L+1)x(2L+1) neighbourhood around point (i,j)
  n <- nrow(field)
  m <- ncol(field)
  i0 <- max(1, i - L)
  i1 <- min(n, i + L)
  j0 <- max(1, j - L)
  j1 <- min(m, j + L)

  neigh <- field[i0:i1, j0:j1]
  if (mode == 'max') {
    return(max(neigh, na.rm = TRUE))
  } else {
    return(min(neigh, na.rm = TRUE))
  }
}

SLX_components <- function(analysis, forecast, L, delta = 0.0) {
  # Identify extreme points (global extremes for simplicity)
  ob_max_pts <- local_extreme_indices(analysis, 'max', delta)
  ob_min_pts <- local_extreme_indices(analysis, 'min', delta)
  fc_max_pts <- local_extreme_indices(forecast, 'max', delta)
  fc_min_pts <- local_extreme_indices(forecast, 'min', delta)
  #browser()
  # Helper to compute average score over a list of points
  avg_score <- function(pts, ob_field, fc_field, mode) {
    if (nrow(pts) == 0) {
      return(NA)
    }

    scores <- c()
    for (idx in 1:nrow(pts)) {
      i <- pts[idx, 1]
      j <- pts[idx, 2]

      if (mode == 'ob_max') {
        ob <- ob_field[i, j]
        phi <- neighbourhood_view(fc_field, i, j, L, 'max')
      } else if (mode == 'ob_min') {
        ob <- ob_field[i, j]
        phi <- neighbourhood_view(fc_field, i, j, L, 'min')
      } else if (mode == 'fc_max') {
        ob <- neighbourhood_view(ob_field, i, j, L, 'max')
        phi <- fc_field[i, j]
      } else if (mode == 'fc_min') {
        ob <- neighbourhood_view(ob_field, i, j, L, 'min')
        phi <- fc_field[i, j]
      } else {
        stop("Invalid mode")
      }

      scores <- c(scores, S_score(ob, phi))
    }
    return(mean(scores, na.rm = TRUE))
  }

  s_ob_max <- avg_score(ob_max_pts, analysis, forecast, 'ob_max')
  s_ob_min <- avg_score(ob_min_pts, analysis, forecast, 'ob_min')
  s_fc_max <- avg_score(fc_max_pts, analysis, forecast, 'fc_max')
  s_fc_min <- avg_score(fc_min_pts, analysis, forecast, 'fc_min')

  combined <- mean(c(s_ob_max, s_ob_min, s_fc_max, s_fc_min), na.rm = TRUE)

  return(list(
    combined = combined,
    s_ob_max = s_ob_max,
    s_ob_min = s_ob_min,
    s_fc_max = s_fc_max,
    s_fc_min = s_fc_min
  ))
}

# ---- Synthetic case (Fig 7) ----
generate_fields <- function(N = 50, block_size = 10, b = 5, d = 15, c = 5) {
  analysis <- matrix(b, nrow = N, ncol = N)
  forecast <- matrix(d, nrow = N, ncol = N)

  # Observed field: two blocks at the top
  # R uses 1-based indexing, so adjust coordinates from Python
  a1_rows <- 16:(15 + block_size)
  a1_cols <- 11:(10 + block_size)
  a2_rows <- 16:(15 + block_size)
  a2_cols <- 31:(30 + block_size)

  analysis[a1_rows, a1_cols] <- b + c
  analysis[a2_rows, a2_cols] <- b + c

  # Forecast field: two blocks at the bottom
  f1_rows <- 31:(30 + block_size)
  f1_cols <- 11:(10 + block_size)
  f2_rows <- 31:(30 + block_size)
  f2_cols <- 31:(30 + block_size)

  forecast[f1_rows, f1_cols] <- d + c
  forecast[f2_rows, f2_cols] <- d + c

  return(list(analysis = analysis, forecast = forecast))
}

plot_fields <- function(analysis, forecast) {
  # Simple plotting function
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

  image(1:ncol(analysis), 1:nrow(analysis), t(analysis), 
        #col = inferno(20),main = "Analysis Field",
        col = heat.colors(20), main = "Analysis Field",
        xlab = "Column", ylab = "Row")

  image(1:ncol(forecast), 1:nrow(forecast), t(forecast), 
        col = heat.colors(20), main = "Forecast Field",
        xlab = "Column", ylab = "Row")

  par(mfrow = c(1, 1))
}

# ---- Main execution ----
cat("=== Simple R SLX Implementation ===\n")

# Test with a single case first
cat("Testing with d=15...\n")
fields <- generate_fields(d = 15)
analysis <- fields$analysis
forecast <- fields$forecast

source("read_data_for_validation_slx.R")
analysis <- obs_field
forecast <- fc_field_5

cat("Analysis field range:", range(analysis), "\n")
cat("Forecast field range:", range(forecast), "\n")

# Plot the fields
# Test SLX for a few L values
cat("\nTesting SLX for different L values:\n")
L_test <- c(0, 5, 10, 15, 20,25)

results <- list()
print("Testing perturbation d=5")
rows <- c()
for (L in L_test) {
  result <- SLX_components(analysis, fc_field_5, L)
  cat(sprintf("L=%d: SLX=%.4f (ob_max=%.3f, ob_min=%.3f, fc_max=%.3f, fc_min=%.3f)\n", 
              L, result$combined, result$s_ob_max, result$s_ob_min, 
              result$s_fc_max, result$s_fc_min))
    rows <- c(rows, result$combined)
}
results[[as.character(5)]] <- rows

print("Testing perturbation d=10")
rows <- c()
for (L in L_test) {
  result <- SLX_components(analysis, fc_field_10, L)
  cat(sprintf("L=%d: SLX=%.4f (ob_max=%.3f, ob_min=%.3f, fc_max=%.3f, fc_min=%.3f)\n", 
              L, result$combined, result$s_ob_max, result$s_ob_min, 
              result$s_fc_max, result$s_fc_min))
    rows <- c(rows, result$combined)
}

results[[as.character(10)]] <- rows
print("Testing perturbation d=20")
rows <- c()
for (L in L_test) {
  result <- SLX_components(analysis, fc_field_20, L)
  cat(sprintf("L=%d: SLX=%.4f (ob_max=%.3f, ob_min=%.3f, fc_max=%.3f, fc_min=%.3f)\n", 
              L, result$combined, result$s_ob_max, result$s_ob_min, 
              result$s_fc_max, result$s_fc_min))
    rows <- c(rows, result$combined)
}
results[[as.character(20)]] <- rows

print("Testing perturbation d=30")
rows <- c()
for (L in L_test) {
  result <- SLX_components(analysis, fc_field_30, L)
  cat(sprintf("L=%d: SLX=%.4f (ob_max=%.3f, ob_min=%.3f, fc_max=%.3f, fc_min=%.3f)\n", 
              L, result$combined, result$s_ob_max, result$s_ob_min, 
              result$s_fc_max, result$s_fc_min))
    rows <- c(rows, result$combined)
}
results[[as.character(30)]] <- rows

png("slx_validation.png", width = 800, height = 600, res = 150)

plot(L_test, results[["5"]], type = "o", col = "red", pch = 19,lwd=2,
     ylim = c(0, 1), xlab = "Neighbourhood width L (grid points)",
     ylab = "Combined SLX", main = "SLX vs L for varying FC perturbation")
lines(L_test, results[["10"]], type = "b", pch = 17, col = "blue", lwd = 2)
lines(L_test, results[["20"]], type = "b", pch = 17, col = "orange", lwd = 2)
lines(L_test, results[["30"]], type = "b", pch = 17, col = "purple", lwd = 2)
legend("topright", legend = c("d=5", "d=10", "d=20", "d=30"),
       col = c("red", "blue", "orange", "purple"),
       pch = c(19, 17, 15, 18), lwd = 2, cex = 0.7)
grid(col = "gray", lty = 2)
dev.off()

# Now test multiple d values like in Python
cat("\n=== Testing multiple d values with the analytic fake data ===\n")
Ds <- c(5, 10, 15, 20, 25,30, 40, 50, 60)
L_values <- seq(0, 40, by = 5)

results <- list()
for (d in Ds) {
  cat(sprintf("Processing d=%d...\n", d))
  fields <- generate_fields(d = d)
  analysis <- fields$analysis
  forecast <- fields$forecast
  #plot_fields(analysis, forecast)
  #browser()

  rows <- c()
  for (L in L_values) {
    result <- SLX_components(analysis, forecast, L)
    rows <- c(rows, result$combined)
  }
  results[[as.character(d)]] <- rows
}

# Plot results similar to Figure 8
plot(L_values, results[["5"]], type = "o", col = "red", pch = 19,lwd=2,
     ylim = c(0, 1), xlab = "Neighbourhood width L (grid points)",
     ylab = "Combined SLX", main = "Synthetic SLX vs neighbourhood width")
#lines(L_values, results[["10"]], type = "b", pch = 17, col = "blue", lwd = 2)
#lines(L_values, results[["20"]], type = "b", pch = 17, col = "orange", lwd = 2)
#lines(L_values, results[["30"]], type = "b", pch = 17, col = "purple", lwd = 2)
#lines(L_values, results[["40"]], type = "b", pch = 17, col = "brown", lwd = 2)
#legend("topright", legend = c("d=5", "d=10", "d=20", "d=30","d=40"),
#       col = c("red", "blue", "orange", "purple","brown"),
#       pch = c(19, 17, 15, 18), lwd = 2, cex = 0.7)
#grid(col = "gray", lty = 2)

colors <- rainbow(length(Ds))
for (i in 1:length(Ds)) {
  d <- Ds[i]
  lines(L_values, results[[as.character(d)]], type = "o", 
        col = colors[i], pch = 19)
}

legend("bottomright", legend = paste("d =", Ds), 
       col = colors, lty = 1, pch = 19, cex = 0.8)
grid()

cat("\n=== Analysis Complete ===\n")
