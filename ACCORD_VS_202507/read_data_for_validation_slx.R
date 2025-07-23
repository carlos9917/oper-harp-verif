#,CORRECTED,SLX Implementation following Sass (2021)

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
fc_file_opts <- netcdf_opts(proj4_var = "RotatedLatLon_Projection", dy = 0.021479999999996835, dx = 0.03607000000000049 ) #, param_find = list("tp"))
cat("Loading forecast data...\n")
precip_fc_5 <- read_grid(file_name="dummy_data/synthetic_tp_forecast_5.nc", parameter="tp_forecast", file_format="netcdf", file_format_opts = fc_file_opts)
precip_fc_10 <- read_grid(file_name="dummy_data/synthetic_tp_forecast_10.nc", parameter="tp_forecast", file_format="netcdf", file_format_opts = fc_file_opts)
precip_fc_20 <- read_grid(file_name="dummy_data/synthetic_tp_forecast_20.nc", parameter="tp_forecast", file_format="netcdf", file_format_opts = fc_file_opts)
precip_fc_30 <- read_grid(file_name="dummy_data/synthetic_tp_forecast_30.nc", parameter="tp_forecast", file_format="netcdf", file_format_opts = fc_file_opts)

# Load observation data
cat("Loading observation data...\n")
precip_ob <- read_grid(file_name="dummy_data/synthetic_tp_analysis.nc", parameter="tp_analysis", file_format="netcdf", file_format_opts = fc_file_opts)

png("fc_field_5_val.png", width = 600, height = 600, res = 120)
plot_field(precip_fc_5)
dev.off()

png("fc_field_10_val.png", width = 600, height = 600, res = 120)
plot_field(precip_fc_10)
dev.off()

png("fc_field_20_val.png", width = 600, height = 600, res = 120)
plot_field(precip_fc_20)
dev.off()

png("fc_field_30_val.png", width = 600, height = 600, res = 120)
plot_field(precip_fc_30)
dev.off()

png("obs_field_val.png", width = 600, height = 600, res = 120)
plot_field(precip_ob)
dev.off()

# Regrid forecast to observation domain
# Convert to arrays and handle NA values
cat("Converting to arrays and handling NA values...\n")
#obs_field <- as.array(precip_ob)
#fc_field <- as.array(precip_fc)

# Get the dimensions first
dims <- dim(precip_ob)
n_rows <- dims[1]  # 122 in your case
n_cols <- dims[2]  # 117 in your case
# Extract as numeric vector and reshape
obs_values <- as.numeric(precip_ob)
obs_field <- matrix(obs_values, nrow = n_rows, ncol = n_cols)

dims <- dim(precip_fc)
n_rows <- dims[1]  # 122 in your case
n_cols <- dims[2]  # 117 in your case
# Extract as numeric vector and reshape
fc_values <- as.numeric(precip_fc_5)
fc_field_5 <- matrix(fc_values, nrow = n_rows, ncol = n_cols)

fc_values <- as.numeric(precip_fc_10)
fc_field_10 <- matrix(fc_values, nrow = n_rows, ncol = n_cols)

fc_values <- as.numeric(precip_fc_20)
fc_field_20 <- matrix(fc_values, nrow = n_rows, ncol = n_cols)

fc_values <- as.numeric(precip_fc_30)
fc_field_30 <- matrix(fc_values, nrow = n_rows, ncol = n_cols)
