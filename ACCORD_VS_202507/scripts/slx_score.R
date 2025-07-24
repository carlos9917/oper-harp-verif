# Simple R implementation following the Python approach
# Based on test_simple_domain.py
# This one gets rid of the issues in the local search
# and works with the simple test case in Bent's paper
library(here)

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

plot_fields <- function(analysis, forecast) {
  # Simple plotting function
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

  image(1:ncol(analysis), 1:nrow(analysis), t(analysis), 
        col = heat.colors(20), main = "Analysis Field",
        xlab = "X", ylab = "Y",useRaster=TRUE)

  image(1:ncol(forecast), 1:nrow(forecast), t(forecast), 
        col = heat.colors(20), main = "Forecast Field",
        xlab = "X", ylab = "Y",useRaster=TRUE)

  par(mfrow = c(1, 1))
}

# ---- Main execution ----
cat("=== Simple R SLX Implementation ===\n")
source(paste0(here::here(),"/ACCORD_VS_202507/scripts/load_field_data_slx.R"))
#source("scripts/load_field_data.R")
analysis <- obs_field
forecast <- fc_field

cat("Analysis field range:", range(analysis), "\n")
cat("Forecast field range:", range(forecast), "\n")

# Plot the fields
# plot_fields(analysis, forecast)

# Test SLX for a few L values
L_test <- c(0, 5, 10, 15, 20,25,30)

cat("\nTesting SLX for different L values:\n")
print(L_test)


# Initialize vectors
L_values <- c()
SLX_values <- c()
ob_max_values <- c()
ob_min_values <- c()
fc_max_values <- c()
fc_min_values <- c()

# Collect results
for (L in L_test) {
  result <- SLX_components(analysis, forecast, L)
  L_values <- c(L_values, L)
  SLX_values <- c(SLX_values, result$combined)
  ob_max_values <- c(ob_max_values, result$s_ob_max)
  ob_min_values <- c(ob_min_values, result$s_ob_min)
  fc_max_values <- c(fc_max_values, result$s_fc_max)
  fc_min_values <- c(fc_min_values, result$s_fc_min)
  cat(sprintf("L=%d: SLX=%.4f (ob_max=%.3f, ob_min=%.3f, fc_max=%.3f, fc_min=%.3f)\n",
              L, result$combined, result$s_ob_max, result$s_ob_min,
              result$s_fc_max, result$s_fc_min))
}

# Plot 1: SLX (Combined) Only
png("slx_overall.png", width = 800, height = 600, res = 150)
plot(L_values, SLX_values, type = "b", pch = 19, col = "darkgreen",lwd = 3,
     xlab = "Neighborhood size (L)", ylab = "SLX (Combined)", main = "Overall SLX Score vs Neighborhood Size",
     ylim = c(0, 1), cex = 1.2)
grid(col = "gray", lty = 2)
# Add performance zones
abline(h = 0.8, col = "green", lty = 2, lwd = 2)
abline(h = 0.6, col = "orange", lty = 2, lwd = 2)
abline(h = 0.4, col = "red", lty = 2, lwd = 2)
text(max(L_values) * 0.7, 0.9, "Excellent", col = "green", cex = 0.8)
text(max(L_values) * 0.7, 0.7, "Good", col = "orange", cex = 0.8)
text(max(L_values) * 0.7, 0.5, "Moderate", col = "red", cex = 0.8)

dev.off()

# Plot 2: Other Scores Together
png("slx_component.png", width = 800, height = 600, res = 150)

plot(L_values, ob_max_values, type = "b", pch = 19, col = "red",lwd = 2,
xlab = "Neighborhood size (L)", ylab = "Component SLX Score",
     main = "SLX Component Scores",
     ylim = c(0, 1))

lines(L_values, ob_min_values, type = "b", pch = 17, col = "blue", lwd = 2)
lines(L_values, fc_max_values, type = "b", pch = 15, col = "orange", lwd = 2)
lines(L_values, fc_min_values, type = "b", pch = 18, col = "purple", lwd = 2)
legend("topright", legend = c("Obs Max", "Obs Min", "FC Max", "FC Min"),
       col = c("red", "blue", "orange", "purple"),
       pch = c(19, 17, 15, 18), lwd = 2, cex = 0.7)
grid(col = "gray", lty = 2)
dev.off()

#for (L in L_test) {
#  result <- SLX_components(analysis, forecast, L)
#  cat(sprintf("L=%d: SLX=%.4f (ob_max=%.3f, ob_min=%.3f, fc_max=%.3f, fc_min=%.3f)\n", 
#              L, result$combined, result$s_ob_max, result$s_ob_min, 
#              result$s_fc_max, result$s_fc_min))
#}
#

cat("\n=== SLX ANALYSIS COMPLETE ===\n")
cat("Key findings:\n")
cat(sprintf("- Peak SLX score: %.3f at L=%d\n", max(SLX_values), L_values[which.max(SLX_values)]))
cat(sprintf("- Score improvement from L=0 to optimal: %.3f\n", max(SLX_values) - SLX_values[1]))

