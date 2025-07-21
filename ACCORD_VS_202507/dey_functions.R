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
window_mean <- function(mat, k) {
  if (k == 0) return(mat)

  ny <- nrow(mat)
  nx <- ncol(mat)
  result <- matrix(NA, ny, nx)

  for (i in 1:ny) {
    for (j in 1:nx) {
      # Define window bounds
      i_min <- max(1, i - k)
      i_max <- min(ny, i + k)
      j_min <- max(1, j - k)
      j_max <- min(nx, j + k)

      # Calculate mean over window
      result[i, j] <- mean(mat[i_min:i_max, j_min:j_max], na.rm = TRUE)
    }
  }

  return(result)
}

# Cumulative sum approach (integral image) - O(N*M) after preprocessing
window_mean_cumsum <- function(mat, k) {
  if (k == 0) return(mat)

  ny <- nrow(mat)
  nx <- ncol(mat)

  # Handle NA values
  mat_clean <- mat
  na_mask <- is.na(mat)
  mat_clean[na_mask] <- 0

  # Create padded matrix for boundary handling
  padded <- matrix(0, ny + 2*k, nx + 2*k)
  padded[(k+1):(ny+k), (k+1):(nx+k)] <- mat_clean

  # Calculate cumulative sum (integral image)
  cumsum_mat <- matrix(0, ny + 2*k + 1, nx + 2*k + 1)

  # Fill cumulative sum matrix
  for (i in 2:(ny + 2*k + 1)) {
    for (j in 2:(nx + 2*k + 1)) {
      cumsum_mat[i, j] <- padded[i-1, j-1] + 
                         cumsum_mat[i-1, j] + 
                         cumsum_mat[i, j-1] - 
                         cumsum_mat[i-1, j-1]
    }
  }

  # Calculate window means using integral image
  result <- matrix(NA, ny, nx)
  window_area <- (2*k + 1)^2

  for (i in 1:ny) {
    for (j in 1:nx) {
      # Calculate indices for integral image lookup
      i1 <- i  # top-left corner in cumsum_mat coordinates
      j1 <- j
      i2 <- i + 2*k + 1  # bottom-right corner
      j2 <- j + 2*k + 1

      # Calculate sum using integral image formula
      window_sum <- cumsum_mat[i2, j2] - 
                   cumsum_mat[i1, j2] - 
                   cumsum_mat[i2, j1] + 
                   cumsum_mat[i1, j1]

      result[i, j] <- window_sum / window_area
    }
  }

  # Restore NA values where original data had NAs
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
    #f1_bar <- window_mean(f1, S)
    #f2_bar <- window_mean(f2, S)
    # optimized
    f1_bar <- window_mean_cumsum(f1, S)
    f2_bar <- window_mean_cumsum(f2, S)

    #f1_bar <- cumsum_2d(f1, S)
    #f2_bar <- cumsum_2d(f2, S)

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
