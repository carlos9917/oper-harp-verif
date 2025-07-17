# Enhanced sass_score_extend.R - MODULAR VERSION WITH SEPARATE CPP FILES
library(harp)
library(harpVis)
library(harpIO)
library(harpSpatial)
library(dplyr)
library(hdf5r)
library(ggplot2)
library(gridExtra)
library(viridis)
library(scico)
library(Rcpp)
library(rlang)

# File paths (adjust as needed)
ob_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/radar_dmi/sqpe/kavrrad_1h/2025/07"
fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/dini_eps/2025070200"


#local DMI data
ob_file_path <- "/ec/res4/hpcperm/nhd/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/radar_dmi/sqpe/kavrrad_1h/2025/07"
fc_file_path <- "/ec/res4/hpcperm/nhd/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/dini_eps/2025070200"

fc_file <- paste0(fc_file_path,"/tp_ekmi_002_mbr001.grib2")
ob_file <-  paste0(ob_file_path,"/202507020200.kavrRAD.01.h5")

#INCA and CLAEF data
ob_file_path <- "/hpcperm/kmek/obs/INCAPlus_1h/inca/2025/07/01"
fc_file_path <- "/hpcperm/kmek/models/CLAEF1k/20250701/00"
ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507011600.nc")
fc_file <- paste0(fc_file_path,"/CLAEF00+0016:00.grb")

veri_time <- "202507011600"

#### LOAD MODULAR C++ FUNCTIONS ####

cat("=== COMPILING C++ FUNCTIONS ===\n")

# Check if C++ files exist
cpp_files <- c("extrema_detection.cpp", "neighborhood_functions.cpp", 
               "score_functions.cpp", "batch_calculations.cpp")

missing_files <- cpp_files[!file.exists(cpp_files)]
if (length(missing_files) > 0) {
  stop(paste("Missing C++ files:", paste(missing_files, collapse=", ")))
}

# Compile each C++ file separately
cat("Compiling extrema detection functions...\n")
sourceCpp("extrema_detection.cpp")

cat("Compiling neighborhood functions...\n")
sourceCpp("neighborhood_functions.cpp")

cat("Compiling score functions...\n")
sourceCpp("score_functions.cpp")

cat("Compiling batch calculation functions...\n")
sourceCpp("batch_calculations.cpp")

cat("All C++ functions compiled successfully!\n")

# Test that functions are available
cat("Testing C++ function availability...\n")
test_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow=3, ncol=3)
test_result <- find_local_extrema_cpp(test_matrix, "max", 0.0)
cat(sprintf("Test successful: found %d extrema in test matrix\n", nrow(test_result)))

cat("=== STARTING DATA LOADING ===\n")
cat("Loading forecast data...\n")
# Read data
#precip_fc <- harpIO::read_grid(paste0(fc_file_path,"/tp_ekmi_002_mbr001.grib2"), "Pcp")
# DINI
#precip_fc <- harpIO::read_grid(fc_file, "Pcp")

# CLAEF
fc_file_opts     <- grib_opts( param_find = setNames( list(list(key = 'indicatorOfParameter', value = 61)), "pcp"))  
precip_fc <- read_grid(file_name=fc_file, parameter="pcp", dttm="20250701", file_format="grib", file_format_opts = fc_file_opts)


dom_fc <- get_domain(precip_fc)
cat("Forecast data loaded successfully\n")

cat("Loading observation data...\n")

############################ DMI radar obs
#this works for DMI radar
#precip_ob <- harpIO::read_grid(
#  ob_file,"Pcp",
#  hdf5_opts = hdf5_opts(data_path = "/pcp/data1/data", odim = FALSE, meta = TRUE)
#)
#############################################

############################ INCA

parameter_inca <- "acc1h"
ob_file_opts     <- netcdf_opts(proj4_var  = "lambert_conformal_conic", param_find = list2(!!parameter_inca := "RR"))
precip_ob <- read_grid(file_name=ob_file, parameter="RR", dttm=veri_time, file_format="netcdf", file_format_opts = ob_file_opts)


#############################################

dom_ob <- get_domain(precip_ob)
cat("Observation data loaded successfully\n")

cat("Regridding forecast to observation domain...\n")
precip_fc_regrid <- geo_regrid(precip_fc, dom_ob)
cat("Regridding completed\n")

# Convert to arrays
cat("Converting to arrays and handling NA values...\n")
obs_field <- as.array(precip_ob)
fc_field <- as.array(precip_fc_regrid)

# Handle NA values
obs_field[is.na(obs_field)] <- 0
fc_field[is.na(fc_field)] <- 0

cat(sprintf("Grid dimensions: %d x %d\n", nrow(obs_field), ncol(obs_field)))
cat("Data preparation completed\n")

#### OPTIMIZED SLX IMPLEMENTATION ####

# Main SLX calculation function using modular C++ functions
calculate_slx_sass_modular <- function(obs, forecast, neighbourhood_sizes=c(0, 1, 3, 5, 7, 9), 
                                      tolerance=0.0, k=0.1, A=4.0) {

  cat("\n=== FINDING LOCAL EXTREMA (MODULAR C++) ===\n")

  # Find local extrema using modular C++ functions
  start_time <- Sys.time()

  cat("Finding observed maxima...\n")
  obs_maxima <- find_local_extrema_cpp(obs, "max", tolerance)

  cat("Finding observed minima...\n")
  obs_minima <- find_local_extrema_cpp(obs, "min", tolerance)

  cat("Finding forecast maxima...\n")
  fc_maxima <- find_local_extrema_cpp(forecast, "max", tolerance)

  cat("Finding forecast minima...\n")
  fc_minima <- find_local_extrema_cpp(forecast, "min", tolerance)

  extrema_time <- Sys.time() - start_time
  cat(sprintf("Extrema detection completed in %.2f seconds\n", as.numeric(extrema_time)))

  cat("\nExtrema detection results:\n")
  cat(sprintf("Observed maxima: %d\n", nrow(obs_maxima)))
  cat(sprintf("Observed minima: %d\n", nrow(obs_minima)))
  cat(sprintf("Forecast maxima: %d\n", nrow(fc_maxima)))
  cat(sprintf("Forecast minima: %d\n", nrow(fc_minima)))

  results <- list()

  cat("\n=== CALCULATING SLX SCORES (MODULAR C++) ===\n")

  for (L in neighbourhood_sizes) {
    cat(sprintf("\nProcessing neighborhood size L = %d\n", L))
    start_L_time <- Sys.time()

    # Calculate component scores using modular C++ functions
    scores_ob_max <- c()
    scores_ob_min <- c()
    scores_fc_max <- c()
    scores_fc_min <- c()

    # SLX_ob_max: How well forecast captures observed maxima
    if (nrow(obs_maxima) > 0) {
      cat(sprintf("  Computing SLX_ob_max for %d observed maxima...\n", nrow(obs_maxima)))
      scores_ob_max <- calculate_component_scores_cpp(obs_maxima, forecast, L, "max", "obs_based", k, A)
    }

    # SLX_ob_min: How well forecast captures observed minima  
    if (nrow(obs_minima) > 0) {
      cat(sprintf("  Computing SLX_ob_min for %d observed minima...\n", nrow(obs_minima)))
      scores_ob_min <- calculate_component_scores_cpp(obs_minima, forecast, L, "min", "obs_based", k, A)
    }

    # SLX_fc_max: How well observed field captures forecast maxima
    if (nrow(fc_maxima) > 0) {
      cat(sprintf("  Computing SLX_fc_max for %d forecast maxima...\n", nrow(fc_maxima)))
      scores_fc_max <- calculate_component_scores_cpp(fc_maxima, obs, L, "max", "fc_based", k, A)
    }

    # SLX_fc_min: How well observed field captures forecast minima
    if (nrow(fc_minima) > 0) {
      cat(sprintf("  Computing SLX_fc_min for %d forecast minima...\n", nrow(fc_minima)))
      scores_fc_min <- calculate_component_scores_cpp(fc_minima, obs, L, "min", "fc_based", k, A)
    }

    # Calculate component averages
    slx_ob_max <- if (length(scores_ob_max) > 0) mean(scores_ob_max, na.rm=TRUE) else 0.0
    slx_ob_min <- if (length(scores_ob_min) > 0) mean(scores_ob_min, na.rm=TRUE) else 0.0
    slx_fc_max <- if (length(scores_fc_max) > 0) mean(scores_fc_max, na.rm=TRUE) else 0.0
    slx_fc_min <- if (length(scores_fc_min) > 0) mean(scores_fc_min, na.rm=TRUE) else 0.0

    # Overall SLX score (Equation 8)
    slx_total <- 0.25 * (slx_ob_max + slx_ob_min + slx_fc_max + slx_fc_min)

    L_time <- Sys.time() - start_L_time
    cat(sprintf("  Completed L=%d: SLX=%.4f (%.2f seconds)\n", L, slx_total, as.numeric(L_time)))


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
      scores_fc_min = scores_fc_min,
      obs_maxima = obs_maxima,
      obs_minima = obs_minima,
      fc_maxima = fc_maxima,
      fc_minima = fc_minima
    )
    #plot_field_with_extrema(L,obs, obs_maxima, obs_minima, 
    #                   "Observed Precipitation with Extrema")
    #plot_field_with_extrema(L,forecast, fc_maxima, fc_minima, 
    #                   "Forecast Precipitation with Extrema")

    #stop("exit here")
  }

  return(results)
}

#### IMPROVED VISUALIZATIONS ####

# Function to create better field plots
plot_field_with_extrema <- function(L,field_data, extrema_max, extrema_min, title, 
                                   max_color="red", min_color="blue") {
print(paste0("plot for: ",L, " ",title))
step <- "extrema"
field <- "observed"
#L <- 1
filename <- sprintf("plot_%s_%s_L%d.png", step, field, L)
png(filename, width = 800, height = 600, res = 150)
  # Create base field plot
  plot_field(field_data, 
             main = title,
             legend_label = "Precipitation (mm/h)",
             col = viridis(20),
             na_colour = "gray95")

  # Add extrema points if they exist
  if (nrow(extrema_max) > 0) {
    points(extrema_max[, "col"], extrema_max[, "row"], 
           col = max_color, pch = 3, cex = 1.5, lwd = 3)
  }

  if (nrow(extrema_min) > 0) {
    points(extrema_min[, "col"], extrema_min[, "row"], 
           col = min_color, pch = 4, cex = 1.5, lwd = 3)
  }

  # Add legend for extrema
  legend("topright", 
         legend = c("Local Max", "Local Min"), 
         col = c(max_color, min_color), 
         pch = c(3, 4), 
         pt.cex = 1.5, 
         pt.lwd = 3,
         bg = "white")
dev.off()
}

# Calculate SLX scores with timing
cat("\n=== STARTING MODULAR SLX CALCULATION ===\n")
total_start_time <- Sys.time()

slx_results <- calculate_slx_sass_modular(obs_field, fc_field)

total_time <- Sys.time() - total_start_time
cat(sprintf("\n=== TOTAL COMPUTATION TIME: %.2f seconds ===\n", as.numeric(total_time)))

cat("\n=== GENERATING RESULTS TABLE ===\n")
# Display results table
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("SLX RESULTS (Modular C++ Implementation)\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat(sprintf("%-3s %-8s %-8s %-8s %-8s %-8s %-12s\n",
           "L", "SLX", "ob_max", "ob_min", "fc_max", "fc_min", "total_extrema"))
cat(paste(rep("-", 80), collapse=""), "\n")

for (L in names(slx_results)) {
  r <- slx_results[[L]]
  total_extrema <- r$n_obs_max + r$n_obs_min + r$n_fc_max + r$n_fc_min
  cat(sprintf("%-3s %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f %-12d\n",
             L, r$SLX, r$SLX_ob_max, r$SLX_ob_min, r$SLX_fc_max, r$SLX_fc_min, total_extrema))
}

cat("\n=== CREATING VISUALIZATIONS ===\n")
#currently not plotting symbols in right places
## VISUALIZATION 1: Original fields with extrema (L=1)
#cat("Creating field plots with extrema...\n")
#par(mfrow = c(1, 2), mar = c(4, 4, 3, 6))
#
#L1_results <- slx_results[["1"]]
#plot_field_with_extrema(precip_ob, L1_results$obs_maxima, L1_results$obs_minima, 
#                       "Observed Precipitation with Extrema")
#plot_field_with_extrema(precip_fc_regrid, L1_results$fc_maxima, L1_results$fc_minima, 
#                       "Forecast Precipitation with Extrema")

# VISUALIZATION 2: SLX component evolution
cat("Creating SLX evolution plots...\n")
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

L_vals <- as.numeric(names(slx_results))
SLX_vals <- sapply(slx_results, function(x) x$SLX)
SLX_ob_max_vals <- sapply(slx_results, function(x) x$SLX_ob_max)
SLX_ob_min_vals <- sapply(slx_results, function(x) x$SLX_ob_min)
SLX_fc_max_vals <- sapply(slx_results, function(x) x$SLX_fc_max)
SLX_fc_min_vals <- sapply(slx_results, function(x) x$SLX_fc_min)

png("slx_overall.png", width = 800, height = 600, res = 150)
# Plot 1: Overall SLX score
plot(L_vals, SLX_vals, type = "b", pch = 19, col = "darkgreen", lwd = 3,
     xlab = "Neighborhood size (L)", ylab = "SLX Score",
     main = "Overall SLX Score vs Neighborhood Size",
     ylim = c(0, 1), cex = 1.2)
grid(col = "gray", lty = 2)

# Add performance zones
abline(h = 0.8, col = "green", lty = 2, lwd = 2)
abline(h = 0.6, col = "orange", lty = 2, lwd = 2)
abline(h = 0.4, col = "red", lty = 2, lwd = 2)
text(max(L_vals) * 0.7, 0.9, "Excellent", col = "green", cex = 0.8)
text(max(L_vals) * 0.7, 0.7, "Good", col = "orange", cex = 0.8)
text(max(L_vals) * 0.7, 0.5, "Moderate", col = "red", cex = 0.8)

dev.off()

# Plot 2: Component scores
png("slx_compnents.png", width = 800, height = 600, res = 150)
plot(L_vals, SLX_ob_max_vals, type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "Neighborhood size (L)", ylab = "Component SLX Score",
     main = "SLX Component Scores",
     ylim = c(0, 1))
lines(L_vals, SLX_ob_min_vals, type = "b", pch = 17, col = "blue", lwd = 2)
lines(L_vals, SLX_fc_max_vals, type = "b", pch = 15, col = "orange", lwd = 2)
lines(L_vals, SLX_fc_min_vals, type = "b", pch = 18, col = "purple", lwd = 2)
legend("topright", legend = c("Obs Max", "Obs Min", "FC Max", "FC Min"),
       col = c("red", "blue", "orange", "purple"),
       pch = c(19, 17, 15, 18), lwd = 2, cex = 0.7)
grid(col = "gray", lty = 2)

dev.off()
# Plot 3: Extrema counts
extrema_counts <- sapply(slx_results, function(x) x$n_obs_max + x$n_obs_min + x$n_fc_max + x$n_fc_min)
plot(L_vals, extrema_counts, type = "b", pch = 19, col = "darkblue", lwd = 2,
     xlab = "Neighborhood size (L)", ylab = "Total Number of Extrema",
     main = "Extrema Count vs Neighborhood Size")
grid(col = "gray", lty = 2)

## Plot 4: Score function visualization
#phi_range <- seq(0, 6, 0.1)
#ob_values <- c(0.05, 1.0, 2.0, 4.0)
#colors <- c("purple", "blue", "green", "red")
#
#plot(phi_range, sapply(phi_range, function(x) score_function_cpp(x, ob_values[1])), 
#     type = "l", col = colors[1], lwd = 2,
#     xlab = "Forecast Value (mm)", ylab = "Score",
#     main = "SLX Score Function", ylim = c(0, 1))
#
#for (i in 2:length(ob_values)) {
#  lines(phi_range, sapply(phi_range, function(x) score_function_cpp(x, ob_values[i])), 
#        col = colors[i], lwd = 2)
#}
#
#abline(v = 0.1, col = "gray", lty = 2, lwd = 1)
#legend("topright", 
#       legend = c(paste("obs =", ob_values), "k = 0.1"),
#       col = c(colors, "gray"), 
#       lty = c(rep(1, 4), 2), 
#       lwd = 2, cex = 0.7)
#grid(col = "gray", lty = 2)
#
## Reset plotting parameters
#par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Key findings:\n")
cat(sprintf("- Peak SLX score: %.3f at L=%d\n", max(SLX_vals), L_vals[which.max(SLX_vals)]))
cat(sprintf("- Score improvement from L=0 to optimal: %.3f\n", max(SLX_vals) - SLX_vals[1]))
cat("- This indicates the spatial tolerance that gives best verification performance\n")
cat(sprintf("- Total computation time: %.2f seconds\n", as.numeric(total_time)))

# Optional: Clean up compiled objects
cat("\nCleaning up temporary compilation files...\n")
temp_files <- list.files(pattern = "^(sourceCpp_|file).*\\.(cpp|o|so|dll)$")
if (length(temp_files) > 0) {
  file.remove(temp_files)
  cat(sprintf("Removed %d temporary files\n", length(temp_files)))
}

cat("Script completed successfully!\n")
