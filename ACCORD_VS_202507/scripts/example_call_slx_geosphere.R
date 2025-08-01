# Example call to the new score
#Key Arguments:
#
#   * scores = "SLX": This tells verify_spatial to execute the SLX score calculation.
#   * window_sizes = c(1, 3, 5, 9, 17, 33): For the SLX score, this vector provides the scales at which to compute
#     the Structure, Amplitude, and Location components.
#   * The file paths, templates, and options are adapted directly from your load_field_data_slx.R script to read the
#     same data.
#   * verif_domain = ob_domain and fc_interp_method = "closest": These arguments ensure the forecast is
#     automatically regridded to the observation grid, just as it was done in the original script.
#   * return_data = TRUE: This makes the function return the score tables as a list.

library(harp)
library(harpIO)
library(harpSpatial)
library(rlang)
parameter_inca <- "accrr1h"

# 1. Set up file format options based on the original script
ob_file_path <- "/hpcperm/kmek/obs/INCAPlus_1h/inca"
ob_file_template <- "{YYYY}/{MM}/{DD}/INCAPlus_1h_RR_ANA_{YYYY}{MM}{DD}{HH}00.nc"

ob_file_format = "netcdf"
parameter <- "accrr12h"
lead_time <- 12
fc_file_opts <- grib_opts(param_find = setNames(list(use_grib_shortName("tp")), parameter))
#ob_file_opts <- netcdf_opts(proj4_var = "lambert_conformal_conic", param_find = list(acc1h = "RR"))
ob_file_opts <- netcdf_opts(proj4_var = "lambert_conformal_conic", param_find = list2(!!parameter := "RR"))

#veri_time <- "202507121200"
##veri_time <- "202507120000"
init_time <- "202507120000" #this is the time when forecast was initialized
#init_time <- "202507121200" #this is the time when forecast was initialized
veri_time <- strftime(strptime(init_time, "%Y%m%d%H") + (lead_time[[1]] * 3600), "%Y%m%d%H00") #this is when the observation is avail 
fc_dttm  <- init_time
#fc_dttm         <- strftime(strptime(veri_time, "%Y%m%d%H") - (lead_time * 3600), "%Y%m%d%H")

ob_file_name <- generate_filenames(file_path     = ob_file_path,
                                   file_date     = veri_time,
                                   file_template = ob_file_template)

# 2. Define the verification domain by reading the observation domain once.
#    verify_spatial will use this to regrid the forecast automatically.
#precip_ob <- read_grid(file_name=ob_file, parameter="RR", dttm=veri_time, file_format="netcdf", file_format_opts = ob_file_opts)


precip_ob <- read_grid(ob_file_name,
                     parameter        = parameter,
                     is_forecast      = TRUE, #FALSE, Vorsicht!
                     dttm             = veri_time,
                     file_format      = ob_file_format,
                     file_format_opts = ob_file_opts)


ob_domain <- get_domain(precip_ob)
print("Read frid from ob")
#ob_domain <- get_domain(
#  generate_filenames(
#    file_path     = "/hpcperm/kmek/obs/INCAPlus_1h/inca/",
#    file_template = "{YYYY}/{MM}/{DD}/INCAPlus_1h_RR_ANA_{YYYY}{MM}{DD}{HH}00.nc",
#    file_date     = "2025071212"
#  ),
#  file_format = "netcdf"
#)

 # 3. Call verify_spatial with the correct arguments
slx_verification <- verify_spatial(
  dttm               =  fc_dttm, # veri_time, #"202507121200",
  parameter          = parameter, #"accrr1h",#"pcp",
  fcst_model         = "CLAEF00",
  lead_time          = lead_time, #12 ,
  scores             = "SLX",
  window_sizes       = c(0, 5, 10, 15, 20,25,30), # c(1, 3, 5, 9, 17, 33), # These are the scales for SLX
  fc_file_path       = "/hpcperm/kmek/models/CLAEF1k",
  fc_file_template   = "{YYYY}{MM}{DD}/00/{det_model}+{LDT4}:00.grb2",
  fc_file_format     = "grib",
  fc_file_opts       = fc_file_opts,
  fc_interp_method   = "nearest", # Method to regrid forecast to obs grid
  ob_file_path       = ob_file_path,
  ob_file_template   = "{YYYY}/{MM}/{DD}/INCAPlus_1h_RR_ANA_{YYYY}{MM}{DD}{HH}00.nc",
  ob_file_format     = "netcdf",
  ob_file_opts       = ob_file_opts,
  ob_accumulation    = "1h",
  verif_domain       = ob_domain,
  return_data        = TRUE # To get the results back in a variable
)

# 4. Print the results
print(slx_verification$SLX)

