
# Figure 4 Validation and Ensemble Extension
# Based on Dey et al. (2016)
# Recreates Figure 4 and demonstrates ensemble functionality

library(fields)  # for image.plot and visualization
library(viridis) # for color palettes
library(here)
library(RColorBrewer)
RdBu_colors <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
# Source the enhanced functions
source(paste0(here::here(),"/ACCORD_VS_202507/scripts/enhanced_dey_functions.R"))

#source("enhanced_dey_functions.R")

cat("=== DEY AGREEMENT SCALES - FIGURE 4 VALIDATION ===\n")

# Set random seed for reproducibility
set.seed(42)

# Parameters from Dey et al. (2016) - Section 4.1
domain_size <- c(193, 242)  # Domain size from paper
n_members <- 12             # 12 member ensemble
rain_area_size <- 50        # L = 50 grid points
rain_area_center <- c(60, 60)  # Lower left corner at (60, 60)
blob_radius <- 8            # Standard rain blob radius
alpha <- 0.5                # Agreement parameter
S_lim <- 80                 # Maximum scale

cat("Parameters:\n")
cat(sprintf("  Domain size: %d x %d\n", domain_size[1], domain_size[2]))
cat(sprintf("  Ensemble members: %d\n", n_members))
cat(sprintf("  Rain area: %d x %d at (%d, %d)\n", rain_area_size, rain_area_size, 
            rain_area_center[1], rain_area_center[2]))
cat(sprintf("  Blob radius: %d grid points\n", blob_radius))
cat(sprintf("  Alpha: %.1f, S_lim: %d\n", alpha, S_lim))

# Create idealized ensemble (Figure 3 setup)
cat("\n=== CREATING IDEALIZED ENSEMBLE ===\n")
ensemble_fields <- create_idealized_ensemble(
  n_members = n_members,
  domain_size = domain_size,
  rain_area_size = rain_area_size,
  rain_area_center = rain_area_center,
  blob_radius = blob_radius,
  rain_value = 1.0
)

png("dey_all_members.png", width = 1600, height = 1200, res = 150)

image(ensemble_fields[[1]], col = c("white", "blue"),
      main = "All Ensemble Members (Contours)", xlab = "Grid X", ylab = "Grid Y", useRaster = TRUE)
for (i in 2:length(ensemble_fields)) {
  contour(ensemble_fields[[i]], add = TRUE, drawlabels = FALSE, col = rgb(0,0,1,0.3))
}
dev.off()

# Create idealized observation (for SA_mo calculation)
cat("\n=== CREATING IDEALIZED OBSERVATION ===\n")
observation_field <- create_idealized_observation(
  domain_size = domain_size,
  obs_area_size = rain_area_size,  # Same size as ensemble area (well-spread case)
  obs_area_center = rain_area_center,  # Same center (well-spread case)
  blob_radius = blob_radius,
  rain_value = 1.0
)

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

png("dey_figure4_validation.png", width = 1600, height = 1200, res = 150)

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

cat("Visualization saved as: dey_figure4_validation.png\n")

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
