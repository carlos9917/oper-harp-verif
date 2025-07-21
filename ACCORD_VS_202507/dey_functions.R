library(harpSpatial)




# CORRECTED Agreement Scales Implementation
# Based on Dey et al. (2016) - Equations 1, 2, 3

# Similarity measure D (Equation 1 from Dey et al. 2016)
similarity_D_corrected <- function(f1_bar, f2_bar) {
  # Handle dimensions properly
  if (is.matrix(f1_bar) || is.array(f1_bar)) {
    f1 <- as.numeric(f1_bar)
    f2 <- as.numeric(f2_bar)
  } else {
    f1 <- f1_bar
    f2 <- f2_bar
  }

  # Initialize D
  D <- numeric(length(f1))

  # Case 1: Both values are zero
  both_zero <- (f1 == 0) & (f2 == 0)
  D[both_zero] <- 1

  # Case 2: At least one value is non-zero
  not_both_zero <- !both_zero
  if (any(not_both_zero)) {
    numerator <- (f1[not_both_zero] - f2[not_both_zero])^2
    denominator <- f1[not_both_zero]^2 + f2[not_both_zero]^2
    D[not_both_zero] <- numerator / denominator
  }

  # Reshape back to original dimensions if needed
  if (is.matrix(f1_bar) || is.array(f1_bar)) {
    D <- array(D, dim = dim(f1_bar))
  }

  return(D)
}

# Simple window mean function
window_mean_cumsum2d <- function(mat, k) {
  if (k == 0L) return(mat)
  
  # Get dimensions
  ny <- nrow(mat)
  nx <- ncol(mat)
  
  # Handle NA values by setting them to 0 for cumsum calculation
  mat_clean <- mat
  na_mask <- is.na(mat)
  mat_clean[na_mask] <- 0
  
  # Calculate cumulative sum using harpSpatial::cumsum_2d
  cumsum_mat <- cumsum_2d(mat_clean)
  
  # Create output matrix
  result <- matrix(NA, ny, nx)
  
  # Calculate window means using integral image
  for (i in 1:ny) {
    for (j in 1:nx) {
      # Define window bounds
      i_min <- max(1, i - k)
      i_max <- min(ny, i + k)
      j_min <- max(1, j - k)
      j_max <- min(nx, j + k)
      
      # Calculate sum using integral image
      # Sum = cumsum[i_max, j_max] - cumsum[i_min-1, j_max] - cumsum[i_max, j_min-1] + cumsum[i_min-1, j_min-1]
      sum_val <- cumsum_mat[i_max, j_max]
      
      if (i_min > 1) {
        sum_val <- sum_val - cumsum_mat[i_min - 1, j_max]
      }
      if (j_min > 1) {
        sum_val <- sum_val - cumsum_mat[i_max, j_min - 1]
      }
      if (i_min > 1 && j_min > 1) {
        sum_val <- sum_val + cumsum_mat[i_min - 1, j_min - 1]
      }
      
      # Calculate number of points in window
      n_points <- (i_max - i_min + 1) * (j_max - j_min + 1)
      
      # Calculate mean
      result[i, j] <- sum_val / n_points
    }
  }
  
  # Restore NA values where original data had NAs in the center point
  result[na_mask] <- NA
  
  return(result)
}

# Agreement scale calculation (corrected)
agreement_scale_map_corrected <- function(f1, f2, alpha = 0.5, S_lim = 80L) {
  # Ensure we work with regular matrices, not geofield objects
  if (inherits(f1, "geofield")) {
    f1 <- as.array(f1)
  }
  if (inherits(f2, "geofield")) {
    f2 <- as.array(f2)
  }

  ny <- nrow(f1)
  nx <- ncol(f1)

  # Initialize agreement scale matrix with maximum scale
  SA <- matrix(S_lim, ny, nx)

  cat(sprintf("Calculating agreement scales for %d x %d grid...\n", ny, nx))

  # Loop over scales from 0 to S_lim
  for (S in 0:S_lim) {
    if (S %% 10 == 0) {
      cat(sprintf("  Processing scale S = %d/%d\n", S, S_lim))
    }

    # Calculate neighborhood means
    #browser()
    #original brute force implementation
    f1_bar <- window_mean_cumsum2d(f1, S)
    f2_bar <- window_mean_cumsum2d(f2, S)


    # Calculate similarity measure D (Equation 1 from Dey et al. 2016)
    D <- similarity_D_corrected(f1_bar, f2_bar)

    # Calculate agreement criterion threshold (Equation 3 from Dey et al. 2016)
    D_crit <- alpha + (1 - alpha) * S / S_lim

    # Find points where agreement is achieved for the first time
    # CORRECTED: Ensure all objects are regular matrices/arrays
    agreement_achieved <- (D <= D_crit) & !is.na(D)
    first_agreement <- agreement_achieved & (SA == S_lim)

    # Update agreement scale for points achieving agreement for first time
    SA[first_agreement] <- S

    # Early termination if all valid points have found their agreement scale
    #remaining_points <- sum(SA == S_lim & !is.na(f1) & !is.na(f2))
    #if (remaining_points == 0) {
    #  cat(sprintf("  All points converged at scale S = %d\n", S))
    #  break
    #}
  }

  cat("Agreement scale calculation completed!\n")
  return(SA)
}
