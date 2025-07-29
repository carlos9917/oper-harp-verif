# plot_slx_from_sqlite.R

# This script reads SLX verification data from a harp-standard SQLite file
# and creates a component plot similar to the one in the original slx_score.R script.

# Ensure the necessary libraries are installed
#if (!requireNamespace("RSQLite", quietly = TRUE)) {
#  install.packages("RSQLite")
#}
#if (!requireNamespace("dplyr", quietly = TRUE)) {
#  install.packages("dplyr")
#}
#if (!requireNamespace("tidyr", quietly = TRUE)) {
#  install.packages("tidyr")
#}

library(RSQLite)
library(dplyr)
library(tidyr)
library(here)

# --- Configuration ---
SQLITE_FILE <- "harp_spatial_scores.sqlite"
OUTPUT_FILE_COMPONENT <- "slx_component_from_sqlite.png"
OUTPUT_FILE_OVERALL   <- "slx_overall_from_sqlite.png"


SQLITE_FILE <- paste0(here::here(),"/ACCORD_VS_202507/harp_spatial_scores.sqlite")
OUTPUT_FILE_COMPONENT <- paste0(here::here(),"/ACCORD_VS_202507/scripts/slx_component_from_sqlite.png")
OUTPUT_FILE_OVERALL <- paste0(here::here(),"/ACCORD_VS_202507/scripts/slx_overall_from_sqlite.png")


FCST_MODEL  <- "CLAEF00" # Specify the model to plot, if more than one is in the db
LEAD_TIME   <- 12        # Specify the lead time (in hours)

# --- Main Script ---

# 1. Connect to the SQLite database
if (!file.exists(SQLITE_FILE)) {
  stop("SQLite file not found: ", SQLITE_FILE)
}
db <- dbConnect(RSQLite::SQLite(), SQLITE_FILE)

# 2. Read the SLX table and filter for the desired case
#    Lead time is stored in seconds in the database.
slx_data <- dbReadTable(db, "SLX") %>%
  filter(
    model == FCST_MODEL,
    leadtime == LEAD_TIME * 3600
  ) %>%
  arrange(scale)

# 3. Disconnect from the database
dbDisconnect(db)

if (nrow(slx_data) == 0) {
  stop("No data found for the specified model and lead time. Please check your configuration.")
}

# 4. Create the component plot
cat("Generating component plot and saving to", OUTPUT_FILE_COMPONENT, "\n")

png(OUTPUT_FILE_COMPONENT, width = 800, height = 600, res = 150)

# Use the same plotting parameters as the original script
plot(slx_data$scale, slx_data$S_OB_MAX, type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "Neighborhood size (L)", ylab = "Component SLX Score",
     main = "SLX Component Scores from SQLite",
     ylim = c(0, 1))

lines(slx_data$scale, slx_data$S_OB_MIN, type = "b", pch = 17, col = "blue", lwd = 2)
lines(slx_data$scale, slx_data$S_FC_MAX, type = "b", pch = 15, col = "orange", lwd = 2)
lines(slx_data$scale, slx_data$S_FC_MIN, type = "b", pch = 18, col = "purple", lwd = 2)

legend("topright", legend = c("Obs Max", "Obs Min", "FC Max", "FC Min"),
       col = c("red", "blue", "orange", "purple"),
       pch = c(19, 17, 15, 18), lwd = 2, cex = 0.7)

grid(col = "gray", lty = 2)

dev.off()

cat("Component plot successfully created.\n\n")

# 5. Create the overall SLX plot
cat("Generating overall plot and saving to", OUTPUT_FILE_OVERALL, "\n")
png(OUTPUT_FILE_OVERALL, width = 800, height = 600, res = 150)

plot(slx_data$scale, slx_data$SLX, type = "b", pch = 19, col = "darkgreen",lwd = 3,
     xlab = "Neighborhood size (L)", ylab = "SLX (Combined)", main = "Overall SLX Score vs Neighborhood Size",
     ylim = c(0, 1), cex = 1.2)
grid(col = "gray", lty = 2)

# Add performance zones
abline(h = 0.8, col = "green", lty = 2, lwd = 2)
abline(h = 0.6, col = "orange", lty = 2, lwd = 2)
abline(h = 0.4, col = "red", lty = 2, lwd = 2)
text(max(slx_data$scale) * 0.7, 0.9, "Excellent", col = "green", cex = 0.8)
text(max(slx_data$scale) * 0.7, 0.7, "Good", col = "orange", cex = 0.8)
text(max(slx_data$scale) * 0.7, 0.5, "Moderate", col = "red", cex = 0.8)

dev.off()
cat("Overall plot successfully created.\n")
