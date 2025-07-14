# Enhanced sass_score_extend.R with visualizations
library(harp)
library(harpVis)
library(harpIO)
library(harpSpatial)
library(dplyr)
library(hdf5r)
library(mmand) # For morphological filtering (local extrema)
library(ggplot2)
library(gridExtra)
library(viridis)
library(scico)
obs_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/radar_dmi/sqpe/kavrrad_1h/2025/07"
#fc_file_path  <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/dini/2025/05"
fc_file_path  <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/dini_eps/2025070200"

# Read model precipitation from GRIB2                  202507020200
precip_fc <- harpIO::read_grid(paste0(fc_file_path,"/tp_ekmi_002_mbr001.grib2"), "Pcp")
dom_fc <-  get_domain(precip_fc)

precip_ob <- harpIO::read_grid(
  paste0(obs_file_path,"/202507020200.kavrRAD.01.h5"),"Pcp",
  hdf5_opts = hdf5_opts(data_path = "/pcp/data1/data", odim = FALSE, meta = TRUE)
)
dom_ob <-  get_domain(precip_ob)

precip_fc_regrid <- geo_regrid(precip_fc, dom_ob)

print("passed regridding fc")

# VISUALIZATION 1: Plot original fields side by side
cat("\n=== VISUALIZATION: Original Fields ===\n")
par(mfrow = c(1, 2))
print("before")
# If you have a geogrid object (from harpIO::read_grid), use it directly:
#plot_field(precip_ob, main = "Observed Precipitation", na_colour = "gray95", col = viridis::viridis(20))
#plot_field(precip_ob,
#  main = "Observed Precipitation",
#  legend_label = "Precipitation (mm/h)",
#  col = viridis(20),
#  na_colour = "gray95")


plot_field(
  precip_ob, 
  palette = scico(256, palette = "cork"), 
  breaks  = seq(0,15,2.5), 
  title = "tot prec over Denmark",
na_colour = "gray95"
)

print("after")
exit("test")
#plot_field(precip_fc_regrid, main = "Forecast Precipitation (Regridded)")
#plot_field(precip_ob, main = "Observed Precipitation", legend_label = "Precipitation (mm/h)")
#plot_field(precip_fc_regrid, main = "Forecast Precipitation (Regridded)", legend_label = "Precipitation (mm/h)")
#plot_field(precip_fc, main = "Forecast Precipitation (Regridded)") #, legend_label = "Precipitation (mm/h)")
plot_field(precip_fc, main = "Forecast Precipitation (Regridded)",
  legend_label = "Precipitation (mm/h)",
  col = viridis(20),
  na_colour = "gray95")




## VISUALIZATION 2: Plot difference field
#par(mfrow = c(1, 1))
## Create difference field manually since plot_field_diff might not be available
#obs_array <- as.array(precip_ob)
#fc_array <- as.array(precip_fc_regrid)
#diff_field <- fc_array - obs_array
#diff_geogrid <- precip_ob  # Copy structure
#diff_geogrid[[1]]$data <- diff_field
#
#print("before plotting diff")
#plot_field(diff_geogrid, main = "Forecast - Observed Difference") #, legend_label = "Difference (mm/h)")

# Convert geogrid objects to simple 2D matrices for calculation
obs_field <- as.array(precip_ob)
fc_field <- as.array(precip_fc_regrid)

# Ensure no NA values, replace with 0
obs_field[is.na(obs_field)] <- 0
fc_field[is.na(fc_field)] <- 0

#### SLX (Sass 2021) Score Calculation ####
# Adapted to use different neighborhood sizes for extrema detection

calculate_slx <- function(obs, forecast, neighbourhood_sizes=c(0, 1, 3, 5, 7, 9), k=0.1, A=4.0, tolerance=0.0) {

  # Helper function to find local extrema with variable neighborhood size
  find_local_extrema <- function(arr, mode='max', tolerance=0.0, L=1) {
    # For L=0, we consider each point as its own extremum
    if (L == 0) {
    indices <- which(!is.na(arr), arr.ind = TRUE)
    if (nrow(indices) == 0) {
    return(matrix(numeric(0), 0, 3, dimnames=list(NULL, c("row", "col", "value"))))
    }
    extrema <- cbind(indices, value = arr[indices])
    return(extrema)
    }
    
    # Create kernel of size (2*L + 1) x (2*L + 1)
    kernel_size <- 2 * L + 1
    kernel <- matrix(1, kernel_size, kernel_size)
    
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
  
  # Helper function for neighbourhood extreme (unchanged)
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
  
  # Helper function for score (unchanged)
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

  results <- list()

  for (L in neighbourhood_sizes) {
    print(paste0("Processing neighborhood size L = ", L))
    
    # Find local extrema using the current neighborhood size L
    obs_maxima <- find_local_extrema(obs, 'max', tolerance, L)
    obs_minima <- find_local_extrema(obs, 'min', tolerance, L)
    fc_maxima <- find_local_extrema(forecast, 'max', tolerance, L)
    fc_minima <- find_local_extrema(forecast, 'min', tolerance, L)
    
    print(paste0("  Found extrema - Obs max: ", nrow(obs_maxima), 
    ", Obs min: ", nrow(obs_minima),
    ", FC max: ", nrow(fc_maxima), 
    ", FC min: ", nrow(fc_minima)))

    # VISUALIZATION 3: Plot extrema for selected neighborhood sizes
    if (L %in% c(0, 1, 3, 7)) {  # Only plot for selected L values to avoid too many plots
      cat(paste0("\\n=== VISUALIZATION: Extrema for L = ", L, " ===\\n"))
      
      # Set up 2x2 plot layout
      par(mfrow = c(2, 2))
      
      # Plot observed maxima
      if (nrow(obs_maxima) > 0) {
        plot_field(precip_ob, main = paste("Observed Maxima (L =", L, ")"))
        points(obs_maxima[, "col"], obs_maxima[, "row"], col = "red", pch = 3, cex = 1.2, lwd = 2)
      } else {
        plot_field(precip_ob, main = paste("Observed Maxima (L =", L, ") - None found"))
      }
      
      # Plot observed minima
      if (nrow(obs_minima) > 0) {
        plot_field(precip_ob, main = paste("Observed Minima (L =", L, ")"))
        points(obs_minima[, "col"], obs_minima[, "row"], col = "blue", pch = 4, cex = 1.2, lwd = 2)
      } else {
        plot_field(precip_ob, main = paste("Observed Minima (L =", L, ") - None found"))
      }
      
      # Plot forecast maxima
      if (nrow(fc_maxima) > 0) {
        plot_field(precip_fc_regrid, main = paste("Forecast Maxima (L =", L, ")"))
        points(fc_maxima[, "col"], fc_maxima[, "row"], col = "red", pch = 3, cex = 1.2, lwd = 2)
      } else {
        plot_field(precip_fc_regrid, main = paste("Forecast Maxima (L =", L, ") - None found"))
      }
      
      # Plot forecast minima
      if (nrow(fc_minima) > 0) {
        plot_field(precip_fc_regrid, main = paste("Forecast Minima (L =", L, ")"))
        points(fc_minima[, "col"], fc_minima[, "row"], col = "blue", pch = 4, cex = 1.2, lwd = 2)
      } else {
        plot_field(precip_fc_regrid, main = paste("Forecast Minima (L =", L, ") - None found"))
      }
    }

    # Calculate scores for observed extrema
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

    # Calculate scores for forecast extrema
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
    n_fc_min = nrow(fc_minima),
    scores_ob_max = scores_ob_max,
    scores_ob_min = scores_ob_min,
    scores_fc_max = scores_fc_max,
    scores_fc_min = scores_fc_min
    )
  }

  return(results)
}

print("Calculating SLX scores with proper neighborhood-based extrema detection")
# Calculate SLX scores
slx_results <- calculate_slx(obs_field, fc_field)

# Display results
cat("\nSLX Results (Sass 2021 Implementation):\n")
cat(paste(rep("=", 70), collapse=""), "\\n")
cat(sprintf("%-3s %-8s %-8s %-8s %-8s %-8s %-12s\\n", 
    "L", "SLX", "ob_max", "ob_min", "fc_max", "fc_min", "extrema_count"))
cat(paste(rep("-", 70), collapse=""), "\\n")

for (L in names(slx_results)) {
    r <- slx_results[[L]]
    total_extrema <- r$n_obs_max + r$n_obs_min + r$n_fc_max + r$n_fc_min
    cat(sprintf("%-3s %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f %-12d\\n",
    L, r$SLX, r$SLX_ob_max, r$SLX_ob_min, r$SLX_fc_max, r$SLX_fc_min, total_extrema))
}

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("Detailed extrema counts by neighborhood size:\n")
cat(sprintf("%-3s %-8s %-8s %-8s %-8s\\n", "L", "obs_max", "obs_min", "fc_max", "fc_min"))
cat(paste(rep("-", 40), collapse=""), "\n")
for (L in names(slx_results)) {
    r <- slx_results[[L]]
    cat(sprintf("%-3s %-8d %-8d %-8d %-8d\n", L, r$n_obs_max, r$n_obs_min, r$n_fc_max, r$n_fc_min))
}

cat("\n")
cat("Note: As neighborhood size L increases, fewer extrema are detected\n")
cat("because larger neighborhoods smooth out local variations.\n")

# VISUALIZATION 4: Summary plots at the end
cat("\n=== FINAL VISUALIZATIONS: Summary Plots ===\n")

# Prepare data for plotting
L_vals <- as.numeric(names(slx_results))
SLX_vals <- sapply(slx_results, function(x) x$SLX)
SLX_ob_max_vals <- sapply(slx_results, function(x) x$SLX_ob_max)
SLX_ob_min_vals <- sapply(slx_results, function(x) x$SLX_ob_min)
SLX_fc_max_vals <- sapply(slx_results, function(x) x$SLX_fc_max)
SLX_fc_min_vals <- sapply(slx_results, function(x) x$SLX_fc_min)
extrema_counts <- sapply(slx_results, function(x) x$n_obs_max + x$n_obs_min + x$n_fc_max + x$n_fc_min)

# Set up 2x2 plot layout for summary
par(mfrow = c(2, 2))

# Plot 1: Overall SLX score vs neighborhood size
plot(L_vals, SLX_vals, type = "b", pch = 19, col = "darkgreen", lwd = 2,
     xlab = "Neighborhood size (L)", ylab = "SLX Score",
     main = "Overall SLX Score vs Neighborhood Size",
     ylim = c(0, 1))
grid()

# Plot 2: Component scores
plot(L_vals, SLX_ob_max_vals, type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "Neighborhood size (L)", ylab = "Component SLX Score",
     main = "SLX Component Scores vs Neighborhood Size",
     ylim = c(0, 1))
lines(L_vals, SLX_ob_min_vals, type = "b", pch = 17, col = "blue", lwd = 2)
lines(L_vals, SLX_fc_max_vals, type = "b", pch = 15, col = "orange", lwd = 2)
lines(L_vals, SLX_fc_min_vals, type = "b", pch = 18, col = "purple", lwd = 2)
legend("topright", legend = c("Obs Max", "Obs Min", "FC Max", "FC Min"),
       col = c("red", "blue", "orange", "purple"), 
       pch = c(19, 17, 15, 18), lwd = 2, cex = 0.8)
grid()

# Plot 3: Number of extrema vs neighborhood size
plot(L_vals, extrema_counts, type = "b", pch = 19, col = "darkblue", lwd = 2,
     xlab = "Neighborhood size (L)", ylab = "Total Number of Extrema",
     main = "Number of Extrema vs Neighborhood Size")
grid()

# Plot 4: Score distribution for L=1 (example)
if ("1" %in% names(slx_results)) {
  all_scores_L1 <- c(slx_results[["1"]]$scores_ob_max, 
                     slx_results[["1"]]$scores_ob_min,
                     slx_results[["1"]]$scores_fc_max, 
                     slx_results[["1"]]$scores_fc_min)
  all_scores_L1 <- all_scores_L1[!is.na(all_scores_L1)]
  
  if (length(all_scores_L1) > 0) {
    hist(all_scores_L1, breaks = 20, col = "lightblue", border = "darkblue",
         xlab = "Individual SLX Scores", ylab = "Frequency",
         main = "Distribution of Individual SLX Scores (L=1)")
    abline(v = mean(all_scores_L1), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste("Mean =", round(mean(all_scores_L1), 3)),
           col = "red", lty = 2, lwd = 2)
  }
} else {
  plot.new()
  text(0.5, 0.5, "No data for L=1", cex = 1.5, adj = 0.5)
}

# Reset plotting parameters
par(mfrow = c(1, 1))

