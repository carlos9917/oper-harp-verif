library(harp)
library(harpIO)
library(harpSpatial)
library(dplyr)
library(hdf5r)

obs_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/radar_dmi/sqpe/kavrrad_1h"

fc_file_path            <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/dini/2025/05"

# Function to calculate Sass score for one threshold
#sass_score_single_threshold <- function(forecast, observation, threshold) {
#  # Binarize forecast and observation at threshold
#  fc_bin <- harp::categorical(forecast, thresholds = threshold)
#  obs_bin <- harp::categorical(observation, thresholds = threshold)
#  
#  # Calculate contingency table components
#  cont <- harp::contingency_table(fc_bin, obs_bin)
#  
#  # Extract hits (H), misses (M), false alarms (F), correct negatives (C)
#  H <- cont$hits
#  M <- cont$misses
#  F <- cont$false_alarms
#  C <- cont$correct_negatives
#  
#  # Sass score formula:
#  # Sass = (H * C - F * M) / ((H + M) * (F + C))
#  numerator <- (H * C) - (F * M)
#  denominator <- (H + M) * (F + C)
#  
#  sass <- ifelse(denominator == 0, NA, numerator / denominator)
#  return(sass)
#}

# Example usage with your data files


# Read model precipitation from GRIB2
precip_fc <- harpIO::read_grid(paste0(fc_file_path,"/tp_2025051300_003.grib2"), "Pcp")
dom_fc <-  get_domain(precip_fc)

# Read radar precipitation from HDF5
#precip_obs <- harpIO::read_grid(
#  paste0(obs_file_path,"/202505130000.kavrRAD.01.h5"),"Pcp",
#  hdf5_opts = hdf5_opts(data_path = "/pcp/data1/data", odim = FALSE, meta = TRUE),
#  transformation = "regrid",
#  transformation_opts = regrid_opts(new_domain=dom_fc)
#)

precip_obs <- harpIO::read_grid(
  paste0(obs_file_path,"/202505130300.kavrRAD.01.h5"),"Pcp",
  hdf5_opts = hdf5_opts(data_path = "/pcp/data1/data", odim = FALSE, meta = TRUE)
)
dom_ob <-  get_domain(precip_ob)

precip_obs_regrid <- geo_regrid(precip_obs, dom_fc)
precip_fc_regrid <- geo_regrid(precip_fc, dom_ob)


print("passed rgridding fc")
# Convert to harp objects
#obs_harp <- harpCore::as_harp_df(precip_obs)
#fc_harp <- harpCore::as_harp_df(precip_fc)

# Define thresholds
thresholds <- c(1, 5, 10, 15, 20, 25)


############calculate contingency tables another way

compute_contingency <- function(forecast, observation, threshold) {
  # Binarize forecast and observation at threshold
  fc_bin <- forecast >= threshold
  obs_bin <- observation >= threshold
  
  # Calculate contingency table elements
  hits <- sum(fc_bin & obs_bin, na.rm = TRUE)
  misses <- sum(!fc_bin & obs_bin, na.rm = TRUE)
  false_alarms <- sum(fc_bin & !obs_bin, na.rm = TRUE)
  correct_negatives <- sum(!fc_bin & !obs_bin, na.rm = TRUE)
  
  return(list(
    hits = hits,
    misses = misses,
    false_alarms = false_alarms,
    correct_negatives = correct_negatives
  ))
}

sass_score_from_contingency <- function(cont) {
  H <- cont$hits
  M <- cont$misses
  F <- cont$false_alarms
  C <- cont$correct_negatives
  
  numerator <- (H * C) - (F * M)
  denominator <- (H + M) * (F + C)
  
  sass <- ifelse(denominator == 0, NA, numerator / denominator)
  return(sass)
}

# Assuming precip_fc and precip_obs are numeric matrices from read_grid()
# Define thresholds
thresholds <- c(1, 5, 10, 15, 20, 25)

# Calculate Sass scores for all thresholds
sass_scores <- sapply(thresholds, function(thr) {
  cont <- compute_contingency(precip_fc_regrid, precip_obs, thr)
  sass_score_from_contingency(cont)
})

# Combine results in a data frame
sass_results <- data.frame(
  Threshold = thresholds,
  Sass_Score = sass_scores
)

print(sass_results)

#######################



## Calculate Sass scores for all thresholds
#sass_scores <- sapply(thresholds, function(thr) {
#  sass_score_single_threshold(fc_harp, obs_harp, thr)
#})
#
## Combine results in a data frame
#sass_results <- data.frame(
#  Threshold = thresholds,
#  Sass_Score = sass_scores
#)
#
#print(sass_results)





#sass_score_single_threshold <- function(forecast, observation, threshold) {
#  # Use harp::verify to get contingency table for threshold
#  verification <- harpSpatial::det_verify(
#    forecast,
#    observation,
#    thresholds = threshold,
#    type = "contingency"
#  )
#  
#  # Extract contingency table values
#  cont <- verification$contingency_table[[1]]
#  
#  # Hits, Misses, False Alarms, Correct Negatives
#  H <- cont$hits
#  M <- cont$misses
#  F <- cont$false_alarms
#  C <- cont$correct_negatives
#  
#  # Sass score formula
#  numerator <- (H * C) - (F * M)
#  denominator <- (H + M) * (F + C)
#  
#  sass <- ifelse(denominator == 0, NA, numerator / denominator)
#  return(sass)
#}
#
## Calculate Sass scores for all thresholds
#sass_scores <- sapply(thresholds, function(thr) {
#  sass_score_single_threshold(fc_harp, obs_harp, thr)
#})
#
## Combine results in a data frame
#sass_results <- data.frame(
#  Threshold = thresholds,
#  Sass_Score = sass_scores
#)
#
#print(sass_results)
