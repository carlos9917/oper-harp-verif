
# Figure 4 Validation and Ensemble Extension
# Based on Dey et al. (2016)
# Recreates Figure 4 and demonstrates ensemble functionality

library(fields)  # for image.plot and visualization
library(viridis) # for color palettes
library(here)
library(RColorBrewer)
library(harp)
library(harpVis)
library(harpIO)
library(harpSpatial)
library(dplyr)
library(hdf5r)
library(ggplot2)
library(gridExtra)
library(scico)
library(rlang)
RdBu_colors <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
# Source the enhanced functions
source(paste0(here::here(),"/ACCORD_VS_202507/scripts/enhanced_dey_functions.R"))


####################### data loading part



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

## data in atos
fc_file_path <- "/hpcperm/kmek/models/CLAEF1k"
ob_file_path <- "/hpcperm/kmek/obs/INCAPlus_1h/inca/2025/07/12"

ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507011600.nc")
fc_file <- paste0(fc_file_path,"/CLAEF00+0016:00.grb")

ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507121200.nc")
fc_file <- paste0(fc_file_path,"/CLAEF00+0012:00.grb2")

fc_file_template <- "{YYYY}{MM}{DD}/00/{fcst_model}{MBR2}+{LDT4}:00.grb2" #for the new files

fcst_model <- "CLAEF"

fc_file_format = "grib"
lead_time <- 12
fc_dttm         <- strftime(strptime(veri_time, "%Y%m%d%H") - (lead_time * 3600), "%Y%m%d%H")
#fc_file_name <- generate_filenames(file_path     = fc_file_path,
#                                   file_date     = fc_dttm,
#                                   lead_time     = lead_time,
#                                   file_template = fc_file_template,
#                                    members       = seq(0,12),         
#                                   fcst_model     = fcst_model)
#

# Load the obs data first, since we will interpolate to this grid

# Load observation data
cat("Loading observation data...\n")
parameter_inca <- "acc1h"
ob_file_opts <- netcdf_opts(proj4_var = "lambert_conformal_conic", param_find = list2(!!parameter_inca := "RR"))
precip_ob <- read_grid(file_name=ob_file, parameter="RR", dttm=veri_time, file_format="netcdf", file_format_opts = ob_file_opts)
dom_ob <- get_domain(precip_ob)
cat("Observation data loaded successfully\n")

# Load forecast data
cat("Loading forecast data...\n")
# for the old files
#fc_file_opts <- grib_opts(param_find = setNames(list(list(key = 'indicatorOfParameter', value = 61)), "pcp"))
# for the new files
fc_file_opts     <- grib_opts( param_find = setNames(list(use_grib_shortName('tp')), "pcp"))

precip_fc <- read_forecast(
                    #fc_file_name,
                    fcst_model          = fcst_model,
                       members        = seq(1, 12),
                     parameter        = "pcp", #"accrr1h",
                     is_forecast      = TRUE,
                     dttm             = fc_dttm,
                     file_format      = fc_file_format,
                     file_template  = fc_file_template,
                     file_format_opts = fc_file_opts,
                     file_path        = fc_file_path,
                     param_defs       = list(),
                     lead_time        = lead_time,
                     transformation      = "regrid",
                     transformation_opts = regrid_opts(new_domain = dom_ob,
                                           method   = "nearest"),
                     return_data   = TRUE 

)

# Get all column names that start with 'CLAEF_mbr'
mbr_cols <- grep("^CLAEF_mbr", names(precip_fc), value = TRUE)

# Extract all geofield matrices into a list
geofields <- lapply(mbr_cols, function(col) {
  precip_fc[[col]][[1]]  # Just [[1]] - this IS the geofield matrix
})

# Now geofields is a list of geofield matrices
names(geofields) <- mbr_cols
cat("Forecast data loaded successfully\n")


#browser()
# Convert to arrays and handle NA values
cat("Converting to arrays and handling NA values...\n")
#obs_field <- as.array(precip_ob)
#fc_field <- as.array(precip_fc_regrid)

# convert the geofields to simple matrices

# OBSERVATION
dims <- dim(precip_ob)
n_rows <- dims[1]
n_cols <- dims[2]
# Extract as numeric vector and reshape
obs_values <- as.numeric(precip_ob)
observation_field <- matrix(obs_values, nrow = n_rows, ncol = n_cols)

# FORECASTS
#dims <- dim(geofields["CLAEF_mbr001"])
#n_rows <- dims[1]
#n_cols <- dims[2]
## Extract as numeric vector and reshape
#fc_values <- as.numeric(precip_fc_regrid)
#fc_field <- matrix(fc_values, nrow = n_rows, ncol = n_cols)


# Get all ensemble member column names
mbr_cols <- grep("^CLAEF_mbr", names(precip_fc), value = TRUE)

# Convert each geofield to a numeric matrix and store in a list
ensemble_fields <- lapply(mbr_cols, function(col) {
  geofield <- precip_fc[[col]][[1]]
  dims <- dim(geofield)
  n_rows <- dims[1]
  n_cols <- dims[2]
  fc_values <- as.numeric(geofield)
  matrix(fc_values, nrow = n_rows, ncol = n_cols)
})

# Name the list elements for clarity
names(ensemble_fields) <- mbr_cols
n_members <- length(mbr_cols)


####################### end data loading part






cat("=== DEY AGREEMENT SCALES  ===\n")


cat("Parameters:\n")
cat(sprintf("  Domain size: %d x %d\n", n_rows, n_cols))
cat(sprintf("  Ensemble members: %d\n", n_members))
cat(sprintf("  Alpha: %.1f, S_lim: %d\n", alpha, S_lim))

# Create idealized ensemble (Figure 3 setup)

png("dey_all_members.png", width = 1600, height = 1200, res = 150)

#image(ensemble_fields[[1]], col = c("white", "blue"),
#      main = "All Ensemble Members (Contours)", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)
#for (i in 2:length(ensemble_fields)) {
#  contour(ensemble_fields[[i]], add = TRUE, drawlabels = FALSE, col = rgb(0,0,1,0.3))
#}

par(mfrow = c(3, 4), mar = c(2, 2, 2, 2))  # Adjust margins as needed

for (i in seq_along(ensemble_fields)) {
  image(
    ensemble_fields[[i]],
    col = c("white", "blue"),
    main = names(ensemble_fields)[i],
    useRaster = TRUE
  )
}
dev.off()

# Calculate SA(mm) - Agreement between ensemble members (Equation 6)
cat("\n=== CALCULATING SA(mm) - ENSEMBLE MEMBER AGREEMENT ===\n")
SA_mm <- calculate_SA_mm(ensemble_fields, alpha = alpha, S_lim = S_lim)

# Calculate SA(mo) - Agreement between members and observation (Equation 7)
cat("\n=== CALCULATING SA(mo) - MEMBER-OBSERVATION AGREEMENT ===\n")
SA_mo <- calculate_SA_mo(ensemble_fields, observation_field, alpha = alpha, S_lim = S_lim)

# Summary statistics
cat("\n=== RESULTS SUMMARY ===\n")
cat(sprintf("SA(mm) - Domain mean: %.1f grid points\n", mean(SA_mm, na.rm = TRUE)))
cat(sprintf("SA(mm) - Min: %.1f, Max: %.1f, SD: %.1f\n", 
            min(SA_mm, na.rm = TRUE), max(SA_mm, na.rm = TRUE), sd(SA_mm, na.rm = TRUE)))

cat(sprintf("SA(mo) - Domain mean: %.1f grid points\n", mean(SA_mo, na.rm = TRUE)))
cat(sprintf("SA(mo) - Min: %.1f, Max: %.1f, SD: %.1f\n", 
            min(SA_mo, na.rm = TRUE), max(SA_mo, na.rm = TRUE), sd(SA_mo, na.rm = TRUE)))

# Create comprehensive visualization
cat("\n=== CREATING VISUALIZATION ===\n")


# Plot the observation in the last slot
image(observation_field, col = c("white", "red"),
      main = "Observation", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)

# Reset plotting layout to default (optional)
par(mfrow = c(1, 1))

png("agreement_members.png", width = 1600, height = 1200, res = 150)

# Set up layout: 3 rows, 3 columns
layout_matrix <- matrix(c(
  1, 2, 3,
  4, 5, 6,
  7, 8, 9
), nrow = 3, byrow = TRUE)
layout(layout_matrix, widths = c(1, 1, 1), heights = c(1, 1, 1))

# Row 1: Show some ensemble members
par(mar = c(3, 3, 3, 1))
image(ensemble_fields[[1]], col = c("white", "blue"), 
      main = "Ensemble Member 1", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)

par(mar = c(3, 3, 3, 1))
image(ensemble_fields[[2]], col = c("white", "blue"), 
      main = "Ensemble Member 2", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)

par(mar = c(3, 3, 3, 1))
image(observation_field, col = c("white", "red"), 
      main = "Observation", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)

# Row 2: Agreement scales
par(mar = c(3, 3, 3, 1))
image(SA_mm, col = plasma(100), 
      main = "SA(mm) - Member Agreement", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)

par(mar = c(3, 3, 3, 1))
image(SA_mo, col = plasma(100), 
      main = "SA(mo) - Member-Obs Agreement", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)

# Difference plot
par(mar = c(3, 3, 3, 1))
SA_diff <- SA_mo - SA_mm
image(SA_diff, col = RdBu_colors,
      main = "SA(mo) - SA(mm)", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)

# Row 3: Analysis plots
# Histogram of SA(mm)
par(mar = c(4, 4, 3, 1))
hist(SA_mm, breaks = 30, col = "lightblue", 
     main = "SA(mm) Distribution", xlab = "Agreement Scale (grid points)", ylab = "Frequency")
abline(v = mean(SA_mm, na.rm = TRUE), col = "red", lwd = 2, lty = 2)

# Histogram of SA(mo)
par(mar = c(4, 4, 3, 1))
hist(SA_mo, breaks = 30, col = "lightgreen", 
     main = "SA(mo) Distribution", xlab = "Agreement Scale (grid points)", ylab = "Frequency")
abline(v = mean(SA_mo, na.rm = TRUE), col = "red", lwd = 2, lty = 2)

# Binned scatter plot (Figure 7 style)
par(mar = c(4, 4, 3, 1))
binned_data <- create_binned_scatter_plot(SA_mm, SA_mo, bin_size = 10, max_scale = S_lim)

plot(binned_data$SA_mm, binned_data$SA_mo, 
     type = "b", pch = 16, col = "blue", lwd = 2,
     xlim = c(0, max(c(binned_data$SA_mm, binned_data$SA_mo), na.rm = TRUE)),
     ylim = c(0, max(c(binned_data$SA_mm, binned_data$SA_mo), na.rm = TRUE)),
     xlab = "SA(mm) - Mean Agreement Scale (grid points)",
     ylab = "SA(mo) - Mean Agreement Scale (grid points)",
     main = "Spread-Skill Relationship")

# Add diagonal line (perfect spread-skill relationship)
abline(0, 1, col = "red", lty = 2, lwd = 2)
legend("topleft", c("Binned Data", "Perfect Spread-Skill"), 
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, pch = c(16, NA))

dev.off()

cat("Visualization saved as: agreement_summary.png\n")

# Create a focused Figure 4 reproduction
cat("\n=== CREATING FIGURE 4 REPRODUCTION ===\n")

# 1. Define the breaks
breaks <- seq(0, 80, by = 10)  # 0, 10, 20, ..., 80

# 2. Create a discrete color palette (e.g., from RColorBrewer)
palette <- brewer.pal(length(breaks) - 1, "RdBu")  # 8 colors

paper_colors <- c(
  "#7B162B", # 0-10
  "#C32F2F", # 10-20
  "#E94B2A", # 20-30
  "#F97B3B", # 30-40
  "#FDBF5B", # 40-50
  "#FEE89A", # 50-60
  "#FEF6C7", # 60-70
  "#FFFFE5"  # 70-80
)

breaks <- seq(0, 80, by = 10)

png("dey_figure4_reproduction.png", width = 800, height = 600, res = 150)

par(mar = c(4, 4, 3, 6))
image(SA_mm,col = paper_colors, breaks = breaks, # col = plasma(100), 
      main = "SA(mm) - Agreement Scales (Figure 4 Reproduction)",
      xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)

# Add colorbar
image.plot(legend.only = TRUE,
           zlim = range(SA_mm, na.rm = TRUE),
           col = paper_colors, #palette, # plasma(100),
           legend.lab = "grid points",
           legend.width = 1.2,
           legend.shrink = 0.8)

dev.off()

cat("Figure 4 reproduction saved as: dey_figure4_reproduction.png\n")

# Validation against paper expectations
cat("\n=== VALIDATION AGAINST PAPER EXPECTATIONS ===\n")
cat("From Dey et al. (2016), Figure 4 description:\n")
cat("- Near the centre of rain area: scales around 10 grid points\n")
cat("- Moving away from precipitation: SA(mm) increases\n")
cat("- Outside rain area: scales represent distance from precipitation\n")

# Check center region (around rain area)
center_region_x <- (rain_area_center[1]):(rain_area_center[1] + rain_area_size)
center_region_y <- (rain_area_center[2]):(rain_area_center[2] + rain_area_size)

# Ensure indices are within bounds
center_region_x <- center_region_x[center_region_x <= domain_size[1]]
center_region_y <- center_region_y[center_region_y <= domain_size[2]]

center_SA_mm <- SA_mm[center_region_x, center_region_y]
center_mean <- mean(center_SA_mm, na.rm = TRUE)

cat(sprintf("\nActual results:\n"))
cat(sprintf("- Center region mean SA(mm): %.1f grid points\n", center_mean))
cat(sprintf("- Domain mean SA(mm): %.1f grid points\n", mean(SA_mm, na.rm = TRUE)))
cat(sprintf("- Expected center ~10, got %.1f - %s\n", center_mean, 
            ifelse(center_mean >= 8 & center_mean <= 15, "GOOD", "CHECK")))

# Check spread-skill relationship for well-spread case
spread_skill_diff <- mean(abs(SA_mo - SA_mm), na.rm = TRUE)
cat(sprintf("- Mean |SA(mo) - SA(mm)|: %.1f (should be small for well-spread)\n", spread_skill_diff))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("1. dey_figure4_validation.png - Comprehensive analysis\n")
cat("2. dey_figure4_reproduction.png - Direct Figure 4 reproduction\n")
cat("3. enhanced_dey_functions.R - Enhanced functions with ensemble support\n")
