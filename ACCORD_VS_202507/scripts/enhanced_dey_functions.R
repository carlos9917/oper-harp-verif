
library(fields)  # for image.plot

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

# Optimized window mean using cumulative sum (integral image)
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
  # Ensure we work with regular matrices
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
    f1_bar <- window_mean_cumsum(f1, S)
    f2_bar <- window_mean_cumsum(f2, S)

    # Calculate similarity measure D (Equation 1 from Dey et al. 2016)
    D <- similarity_D_corrected(f1_bar, f2_bar)

    # Calculate agreement criterion threshold (Equation 3 from Dey et al. 2016)
    D_crit <- alpha + (1 - alpha) * S / S_lim

    # Find points where agreement is achieved for the first time
    agreement_achieved <- (D <= D_crit) & !is.na(D)
    first_agreement <- agreement_achieved & (SA == S_lim)

    # Update agreement scale for points achieving agreement for first time
    SA[first_agreement] <- S
  }

  cat("Agreement scale calculation completed!\n")
  return(SA)
}

# NEW: Ensemble agreement scales (SA_mm) - Equation 6
calculate_SA_mm <- function(ensemble_fields, alpha = 0.5, S_lim = 80L) {
  n_members <- length(ensemble_fields)

  if (n_members < 2) {
    stop("Need at least 2 ensemble members to calculate SA(mm)")
  }

  # Calculate number of member pairs (Np in Equation 5)
  n_pairs <- n_members * (n_members - 1) / 2
  cat(sprintf("Calculating SA(mm) for %d members (%d pairs)...\n", n_members, n_pairs))

  # Store all pairwise agreement scales
  agreement_scales <- list()
  pair_count <- 0

  # Calculate agreement scales for all member pairs
  for (i in 1:(n_members - 1)) {
    for (j in (i + 1):n_members) {
      pair_count <- pair_count + 1
      cat(sprintf("Processing member pair %d-%d (%d/%d)\n", i, j, pair_count, n_pairs))

      scale_map <- agreement_scale_map_corrected(
        ensemble_fields[[i]], 
        ensemble_fields[[j]], 
        alpha = alpha, 
        S_lim = S_lim
      )

      agreement_scales[[pair_count]] <- scale_map
    }
  }

  # Average over all pairs (Equation 6 from Dey et al. 2016)
  SA_mm <- Reduce("+", agreement_scales) / length(agreement_scales)

  cat(sprintf("SA(mm) calculation completed. Mean value: %.2f grid points\n", 
              mean(SA_mm, na.rm = TRUE)))

  return(SA_mm)
}

# NEW: Member-observation agreement scales (SA_mo) - Equation 7
calculate_SA_mo <- function(ensemble_fields, observations, alpha = 0.5, S_lim = 80L) {
  n_members <- length(ensemble_fields)

  cat(sprintf("Calculating SA(mo) for %d members vs observations...\n", n_members))

  # Store agreement scales for all member-observation pairs
  agreement_scales <- list()

  # Calculate agreement scales for each member vs observations
  for (i in 1:n_members) {
    cat(sprintf("Processing member %d/%d vs observations\n", i, n_members))

    scale_map <- agreement_scale_map_corrected(
      ensemble_fields[[i]], 
      observations, 
      alpha = alpha, 
      S_lim = S_lim
    )

    agreement_scales[[i]] <- scale_map
  }

  # Average over all members (Equation 7 from Dey et al. 2016)
  SA_mo <- Reduce("+", agreement_scales) / length(agreement_scales)

  cat(sprintf("SA(mo) calculation completed. Mean value: %.2f grid points\n", 
              mean(SA_mo, na.rm = TRUE)))

  return(SA_mo)
}

# NEW: Create idealized ensemble (Figure 4 setup)
create_idealized_ensemble <- function(n_members = 12, domain_size = c(193, 242), 
                                    rain_area_size = 50, rain_area_center = c(60, 60),
                                    blob_radius = 8, rain_value = 1.0) {

  cat("Creating idealized ensemble...\n")
  cat(sprintf("Domain: %d x %d, Members: %d, Rain area: %d x %d\n", 
              domain_size[1], domain_size[2], n_members, rain_area_size, rain_area_size))

  ensemble_fields <- list()

  # Create each ensemble member
  for (i in 1:n_members) {
    # Initialize domain with zeros
    field <- matrix(0, domain_size[1], domain_size[2])

    # Randomly position rain blob within rain area
    x_center <- rain_area_center[1] + runif(1, 0, rain_area_size)
    y_center <- rain_area_center[2] + runif(1, 0, rain_area_size)

    # Create circular rain blob
    for (x in 1:domain_size[1]) {
      for (y in 1:domain_size[2]) {
        distance <- sqrt((x - x_center)^2 + (y - y_center)^2)
        if (distance <= blob_radius) {
          field[x, y] <- rain_value
        }
      }
    }

    ensemble_fields[[i]] <- field
    cat(sprintf("Member %d: blob center at (%.1f, %.1f)\n", i, x_center, y_center))
  }

  cat("Idealized ensemble created successfully!\n")
  return(ensemble_fields)
}

# NEW: Create idealized observation
create_idealized_observation <- function(domain_size = c(193, 242), 
                                       obs_area_size = 50, obs_area_center = c(60, 60),
                                       blob_radius = 8, rain_value = 1.0) {

  # Initialize domain with zeros
  field <- matrix(0, domain_size[1], domain_size[2])

  # Position observation blob within observation area
  x_center <- obs_area_center[1] + runif(1, 0, obs_area_size)
  y_center <- obs_area_center[2] + runif(1, 0, obs_area_size)

  # Create circular rain blob
  for (x in 1:domain_size[1]) {
    for (y in 1:domain_size[2]) {
      distance <- sqrt((x - x_center)^2 + (y - y_center)^2)
      if (distance <= blob_radius) {
        field[x, y] <- rain_value
      }
    }
  }

  cat(sprintf("Idealized observation: blob center at (%.1f, %.1f)\n", x_center, y_center))
  return(field)
}

# NEW: Binned scatter plot for spread-skill analysis (Figure 7 style)
create_binned_scatter_plot <- function(SA_mm, SA_mo, bin_size = 10, max_scale = 80) {

  # Create bins
  bins <- seq(0, max_scale, by = bin_size)
  n_bins <- length(bins) - 1

  bin_centers <- numeric(n_bins)
  SA_mm_binned <- numeric(n_bins)
  SA_mo_binned <- numeric(n_bins)

  for (i in 1:n_bins) {
    # Find points in this bin
    in_bin <- (SA_mm >= bins[i]) & (SA_mm < bins[i+1]) & !is.na(SA_mm) & !is.na(SA_mo)

    if (sum(in_bin) > 0) {
      bin_centers[i] <- mean(SA_mm[in_bin], na.rm = TRUE)
      SA_mm_binned[i] <- mean(SA_mm[in_bin], na.rm = TRUE)
      SA_mo_binned[i] <- mean(SA_mo[in_bin], na.rm = TRUE)
    } else {
      bin_centers[i] <- NA
      SA_mm_binned[i] <- NA
      SA_mo_binned[i] <- NA
    }
  }

  # Remove empty bins
  valid_bins <- !is.na(bin_centers)

  return(list(
    SA_mm = SA_mm_binned[valid_bins],
    SA_mo = SA_mo_binned[valid_bins],
    bin_centers = bin_centers[valid_bins]
  ))
}
