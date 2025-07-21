
# CORRECTED SLX Implementation following Sass (2021)

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
library(rlang)

#### DATA LOADING SECTION ####
cat("=== STARTING DATA LOADING ===\n")

# File paths and configuration
veri_time <- "202507011600"
veri_time <- "202507121200"
ob_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/INCAPlus_1h/inca/2025/07/01"
ob_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/INCAPlus_1h/inca/2025/07/12"

fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k/20250701/00"
fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k/20250712/00"
fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k"

ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507011600.nc")
fc_file <- paste0(fc_file_path,"/CLAEF00+0016:00.grb")

ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507121200.nc")
fc_file <- paste0(fc_file_path,"/CLAEF00+0012:00.grb2")

fc_file_template <- "{YYYY}{MM}{DD}/00/{det_model}+{LDT4}:00.grb2" #for the new files

fcst_model <- "CLAEF00"

fc_file_format = "grib"
lead_time <- 12
fc_dttm         <- strftime(strptime(veri_time, "%Y%m%d%H") - (lead_time * 3600), "%Y%m%d%H")
fc_file_name <- generate_filenames(file_path     = fc_file_path,
                                   file_date     = fc_dttm,
                                   lead_time     = lead_time,
                                   file_template = fc_file_template,
                                   det_model     = fcst_model)


# Load forecast data
cat("Loading forecast data...\n")
# for the old files
#fc_file_opts <- grib_opts(param_find = setNames(list(list(key = 'indicatorOfParameter', value = 61)), "pcp"))
# for the new files
fc_file_opts     <- grib_opts( param_find = setNames(list(use_grib_shortName('tp')), "pcp"))

#precip_fc <- read_grid(file_name=fc_file, parameter="pcp", dttm="20250701", file_format="grib", file_format_opts = fc_file_opts)
#precip_fc <- read_grid(file_name=fc_file, parameter="pcp", dttm=veri_time, file_format="grib", file_format_opts = fc_file_opts)
#ignoring time to check only
#precip_fc <- read_grid(file_name=fc_file, parameter="pcp", file_format="grib", file_format_opts = fc_file_opts)


#parameter <- "accrr1h"

precip_fc <- read_grid(fc_file_name,
                     parameter        = "pcp", #"accrr1h",
                     is_forecast      = TRUE,
                     dttm             = fc_dttm,
                     file_format      = fc_file_format,
                     file_format_opts = fc_file_opts,
                     param_defs       = list(),
                     lead_time        = lead_time)


dom_fc <- get_domain(precip_fc)
cat("Forecast data loaded successfully\n")

# Load observation data
cat("Loading observation data...\n")
parameter_inca <- "acc1h"
ob_file_opts <- netcdf_opts(proj4_var = "lambert_conformal_conic", param_find = list2(!!parameter_inca := "RR"))
precip_ob <- read_grid(file_name=ob_file, parameter="RR", dttm=veri_time, file_format="netcdf", file_format_opts = ob_file_opts)
dom_ob <- get_domain(precip_ob)
cat("Observation data loaded successfully\n")

# Regrid forecast to observation domain
cat("Regridding forecast to observation domain...\n")
precip_fc_regrid <- geo_regrid(precip_fc, dom_ob)
cat("Regridding completed\n")

# Convert to arrays and handle NA values
cat("Converting to arrays and handling NA values...\n")
obs_field <- as.array(precip_ob)
fc_field <- as.array(precip_fc_regrid)
browser()

#png("obs_field.png", width = 800, height = 600, res = 150)
#plot_field(obs_field)
#dev.off()
#png("fc_field.png", width = 800, height = 600, res = 150)
#plot_field(fc_field)
#dev.off()

# Handle NA values
#obs_field[is.na(obs_field)] <- 0
#fc_field[is.na(fc_field)] <- 0

cat(sprintf("Grid dimensions: %d x %d\n", nrow(obs_field), ncol(obs_field)))
cat("Data preparation completed\n")


source("dey_functions.R")
############################################################
## AGREEMENT-SCALES (Dey et al., 2016) – single pair      ##
############################################################
library(harpSpatial)   # dev branch – optional, see comment

#--- Parameters ----------------------------------------------------------------



#--- Ensemble agreement scales (SA_mm) ----------------------------------------
calculate_SA_mm <- function(ensemble_fields, alpha = 0.5, S_lim = 80L) {
  n_members <- length(ensemble_fields)
  
  if (n_members < 2) {
    stop("Need at least 2 ensemble members to calculate SA(mm)")
  }
  
  # Calculate number of member pairs
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

#--- Member-observation agreement scales (SA_mo) ------------------------------
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

#--- Main execution for single forecast vs observation ------------------------
cat("\n=== CALCULATING AGREEMENT SCALES (Dey et al. 2016) ===\n")

# For single forecast vs observation (SA_mo equivalent)
cat("Calculating agreement scales between forecast and observation...\n")
SA_fo <- agreement_scale_map_corrected(fc_field, obs_field, alpha = alpha, S_lim = S_lim)
print(SA_fo)

browser()
# Summary statistics
cat("\n=== RESULTS SUMMARY ===\n")
cat(sprintf("Domain-mean agreement scale: %.1f grid points\n", mean(SA_fo, na.rm = TRUE)))
cat(sprintf("Minimum agreement scale: %.1f grid points\n", min(SA_fo, na.rm = TRUE)))
cat(sprintf("Maximum agreement scale: %.1f grid points\n", max(SA_fo, na.rm = TRUE)))
cat(sprintf("Standard deviation: %.1f grid points\n", sd(SA_fo, na.rm = TRUE)))

# Create visualization
png("agreement_scales_dey2016.png", width = 1000, height = 800, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 6))

# Plot 1: Observations
image(t(obs_field[nrow(obs_field):1, ]), 
      col = viridis::viridis(100), 
      main = "Observations", 
      xlab = "Grid X", ylab = "Grid Y")

# Plot 2: Forecast
image(t(fc_field[nrow(fc_field):1, ]), 
      col = viridis::viridis(100), 
      main = "Forecast", 
      xlab = "Grid X", ylab = "Grid Y")

# Plot 3: Agreement scales
image(t(SA_fo[nrow(SA_fo):1, ]), 
      col = viridis::plasma(100), 
      main = "Agreement Scales SA(fo)", 
      xlab = "Grid X", ylab = "Grid Y")

# Plot 4: Histogram of agreement scales
hist(SA_fo, breaks = 50, col = "lightblue", 
     main = "Distribution of Agreement Scales",
     xlab = "Agreement Scale (grid points)", 
     ylab = "Frequency")
abline(v = mean(SA_fo, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
legend("topright", "Mean", col = "red", lty = 2, lwd = 2)

dev.off()

cat("\nVisualization saved as 'agreement_scales_dey2016.png'\n")

#--- Optional: For ensemble analysis (if ensemble data available) -------------
# Uncomment and modify if you have ensemble data:
#
# # Convert single forecast to list format for ensemble functions
# ensemble_fields <- list(fc_field)  # Add more members if available
# 
# # Calculate SA(mm) if you have multiple ensemble members
# if (length(ensemble_fields) > 1) {
#   SA_mm <- calculate_SA_mm(ensemble_fields, alpha = alpha, S_lim = S_lim)
#   
#   # Calculate SA(mo)
#   SA_mo <- calculate_SA_mo(ensemble_fields, obs_field, alpha = alpha, S_lim = S_lim)
#   
#   # Spatial spread-skill analysis (as in Dey et al. 2016, Figure 13)
#   # ... add binned scatter plot analysis here
# }

cat("\n=== AGREEMENT SCALES ANALYSIS COMPLETE ===\n")
