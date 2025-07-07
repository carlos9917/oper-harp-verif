library(harp)
library(harpIO)
library(harpSpatial)
library(dplyr)
library(hdf5r)
library(mmand) # For morphological filtering (local extrema)

obs_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/radar_dmi/sqpe/kavrrad_1h"
fc_file_path  <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/dini/2025/05"

# Read model precipitation from GRIB2
precip_fc <- harpIO::read_grid(paste0(fc_file_path,"/tp_2025051300_006.grib2"), "Pcp")
dom_fc <-  get_domain(precip_fc)

precip_ob <- harpIO::read_grid(
  paste0(obs_file_path,"/202505130600.kavrRAD.01.h5"),"Pcp",
  hdf5_opts = hdf5_opts(data_path = "/pcp/data1/data", odim = FALSE, meta = TRUE)
)
dom_ob <-  get_domain(precip_ob)

precip_fc_regrid <- geo_regrid(precip_fc, dom_ob)

print("passed rgridding fc")

# Convert geogrid objects to simple 2D matrices for calculation
# Assuming the structure is [lon, lat, leadtime, member]
obs_field <- as.array(precip_ob)[,,1,1]
fc_field <- as.array(precip_fc_regrid)[,,1,1]

# Ensure no NA values, replace with 0
obs_field[is.na(obs_field)] <- 0
fc_field[is.na(fc_field)] <- 0


############ SLX (Sass 2021) Score Calculation ############

# Based on python implementation from:
# https://github.com/carlos9917/oper-harp-verif/blob/master/ACCORD_VS_202507/docs/review_sass/slx_presentation.qmd

calculate_slx <- function(obs, forecast, neighbourhood_sizes=c(0, 1, 3, 5, 7, 9), k=0.1, A=4.0, tolerance=0.0) {

  # Helper function to find local extrema
  find_local_extrema <- function(arr, mode='max', tolerance=0.0) {
    kernel <- matrix(1, 3, 3)
    if (mode == 'max') {
      filtered <- mmand::dilate(arr, kernel)
      mask <- (arr >= filtered - tolerance) & (arr == filtered)
    } else {
      filtered <- mmand::erode(arr, kernel)
      mask <- (arr <= filtered + tolerance) & (arr == filtered)
    }
    
    indices <- which(mask, arr.ind = TRUE)
    
    if (nrow(indices) == 0) {
      return(matrix(numeric(0), 0, 3, dimnames=list(NULL, c("row", "col", "value"))))
    }
    
    extrema <- cbind(indices, value = arr[indices])
    return(extrema)
  }
  
  # Helper function for neighbourhood extreme
  get_neighbourhood_extreme <- function(arr, i, j, L, mode='max') {
    i_min <- max(1, i - L)
    i_max <- min(nrow(arr), i + L)
    j_min <- max(1, j - L)
    j_max <- min(ncol(arr), j + L)
    
    neighbourhood <- arr[i_min:i_max, j_min:j_max]
    
    if (mode == 'max') {
      return(max(neighbourhood, na.rm = TRUE))
    } else {
      return(min(neighbourhood, na.rm = TRUE))
    }
  }
  
  # Helper function for score
  score_function <- function(phi, ob, k=0.1, A=4.0) {
    if (is.na(phi) || is.na(ob)) return(NA)
    
    if (ob <= k) {
      if (phi <= k) return(1.0)
      else return(max(0.0, 1 - (phi - k) / (A * k)))
    } else {
      if (phi < ob - k) return(phi / (ob - k))
      else if (phi <= ob) return(1.0)
      else return(max(0.0, 1 - (phi - ob) / (A * ob)))
    }
  }

  # Find local extrema
  obs_maxima <- find_local_extrema(obs, 'max', tolerance)
  obs_minima <- find_local_extrema(obs, 'min', tolerance)
  fc_maxima <- find_local_extrema(forecast, 'max', tolerance)
  fc_minima <- find_local_extrema(forecast, 'min', tolerance)

  results <- list()

  for (L in neighbourhood_sizes) {
    scores_ob_max <- if (nrow(obs_maxima) > 0) {
      apply(obs_maxima, 1, function(extr) {
        fc_neighbourhood_max <- get_neighbourhood_extreme(forecast, extr['row'], extr['col'], L, 'max')
        score_function(fc_neighbourhood_max, extr['value'], k, A)
      })
    } else { c() }

    scores_ob_min <- if (nrow(obs_minima) > 0) {
      apply(obs_minima, 1, function(extr) {
        fc_neighbourhood_min <- get_neighbourhood_extreme(forecast, extr['row'], extr['col'], L, 'min')
        score_function(fc_neighbourhood_min, extr['value'], k, A)
      })
    } else { c() }

    scores_fc_max <- if (nrow(fc_maxima) > 0) {
      apply(fc_maxima, 1, function(extr) {
        obs_neighbourhood_max <- get_neighbourhood_extreme(obs, extr['row'], extr['col'], L, 'max')
        score_function(extr['value'], obs_neighbourhood_max, k, A)
      })
    } else { c() }

    scores_fc_min <- if (nrow(fc_minima) > 0) {
      apply(fc_minima, 1, function(extr) {
        obs_neighbourhood_min <- get_neighbourhood_extreme(obs, extr['row'], extr['col'], L, 'min')
        score_function(extr['value'], obs_neighbourhood_min, k, A)
      })
    } else { c() }

    # Calculate component scores
    slx_ob_max <- if (length(scores_ob_max) > 0) mean(scores_ob_max, na.rm=TRUE) else 0.0
    slx_ob_min <- if (length(scores_ob_min) > 0) mean(scores_ob_min, na.rm=TRUE) else 0.0
    slx_fc_max <- if (length(scores_fc_max) > 0) mean(scores_fc_max, na.rm=TRUE) else 0.0
    slx_fc_min <- if (length(scores_fc_min) > 0) mean(scores_fc_min, na.rm=TRUE) else 0.0

    # Overall SLX score
    slx_total <- 0.25 * (slx_ob_max + slx_ob_min + slx_fc_max + slx_fc_min)

    results[[as.character(L)]] <- list(
      SLX = slx_total,
      SLX_ob_max = slx_ob_max,
      SLX_ob_min = slx_ob_min,
      SLX_fc_max = slx_fc_max,
      SLX_fc_min = slx_fc_min,
      n_obs_max = nrow(obs_maxima),
      n_obs_min = nrow(obs_minima),
      n_fc_max = nrow(fc_maxima),
      n_fc_min = nrow(fc_minima)
    )
  }

  return(results)
}

# Calculate SLX scores
slx_results <- calculate_slx(obs_field, fc_field)

# Display results
cat("SLX Results:\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat(sprintf("%-3s %-6s %-7s %-7s %-7s %-7s\n", "L", "SLX", "ob_max", "ob_min", "fc_max", "fc_min"))
cat(paste(rep("-", 60), collapse=""), "\n")

for (L in names(slx_results)) {
    r <- slx_results[[L]]
    cat(sprintf("%-3s %-6.3f %-7.3f %-7.3f %-7.3f %-7.3f\n",
          L, r$SLX, r$SLX_ob_max, r$SLX_ob_min, r$SLX_fc_max, r$SLX_fc_min))
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("Extrema counts:\n")
r_counts <- slx_results[['0']]
cat(sprintf("Observed maxima: %d\n", r_counts$n_obs_max))
cat(sprintf("Observed minima: %d\n", r_counts$n_obs_min))
cat(sprintf("Forecast maxima: %d\n", r_counts$n_fc_max))
cat(sprintf("Forecast minima: %d\n", r_counts$n_fc_min))