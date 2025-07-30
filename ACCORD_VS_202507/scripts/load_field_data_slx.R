library(harp)
library(harpIO)
library(harpSpatial)
library(hdf5r)
library(rlang)
#### DATA LOADING SECTION ####
cat("=== STARTING DATA LOADING ===\n")


# File paths and configuration
veri_time <- "202507011600"
veri_time <- "202507121200"
ob_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/INCAPlus_1h/inca/2025/07/01"
ob_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/INCAPlus_1h/inca/2025/07/12"

#local paths in my laptop
fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k/20250701/00"
fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k/20250712/00"
fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k"

## paths in atos
fc_file_path <- "/hpcperm/kmek/models/CLAEF1k"
ob_file_path <- "/hpcperm/kmek/obs/INCAPlus_1h/inca/2025/07/12"

#ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507011600.nc")
fc_file <- paste0(fc_file_path,"/CLAEF00+0016:00.grb")

ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507121200.nc") #note hardcoded date here, to change
#fc_file <- paste0(fc_file_path,"/CLAEF00+0012:00.grb2") # using fc_file_name below

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

print(fc_file_name)

# Load forecast data
cat("Loading forecast data...\n")
# for the old files
#fc_file_opts <- grib_opts(param_find = setNames(list(list(key = 'indicatorOfParameter', value = 61)), "pcp"))
# for the new files
fc_file_opts     <- grib_opts( param_find = setNames(list(use_grib_shortName('tp')), "pcp"))
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

# Convert to matrices for cleaer handling in the algorithm
cat("Converting to arrays and handling NA values...\n")

# quick and dirty plots
print("plotting obs")
png("obs_field.png", width = 800, height = 600, res = 150)
plot_field(precip_ob,title="OB")
dev.off()

print("plotting fc")
png("fc_field.png", width = 800, height = 600, res = 150)
plot_field(precip_fc_regrid,title="FC")
dev.off() 


# Get the dimensions first
dims <- dim(precip_ob)
n_rows <- dims[1]  
n_cols <- dims[2]  
# Extract as numeric vector and reshape
obs_values <- as.numeric(precip_ob)
obs_field <- matrix(obs_values, nrow = n_rows, ncol = n_cols)

dims <- dim(precip_fc_regrid)
n_rows <- dims[1] 
n_cols <- dims[2] 
# Extract as numeric vector and reshape
fc_values <- as.numeric(precip_fc_regrid)
fc_field <- matrix(fc_values, nrow = n_rows, ncol = n_cols)

