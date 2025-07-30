# Example call to the SLX score using DMI data
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

# 1. Set up file format options based on the original script
#ob_file_path <- "/hpcperm/kmek/obs/INCAPlus_1h/inca"
#ob_file_template <- "{YYYY}/{MM}/{DD}/INCAPlus_1h_RR_ANA_{YYYY}{MM}{DD}{HH}00.nc"

ob_file_path <- "/ec/res4/hpcperm/nhd/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/radar_dmi/sqpe/kavrrad_3h/2025/07"
#ob_file_template <- "{YYYY}/{MM}/{DD}/INCAPlus_1h_RR_ANA_{YYYY}{MM}{DD}{HH}00.nc"
#ob_file_template <- "{YYYY}{MM}{DD}{HH}00.kavrRAD.03.h5"

pcp_accum_period = "3h"  # "1h", "3h", "6h"

parameter            = paste0("Accpcp", pcp_accum_period)
ob_accumulation  <- pcp_accum_period


ob_file_template = switch(
                          pcp_accum_period,
                          "1h"  = "/{YYYY}{MM}{DD}{HH}00.kavrRAD.01.h5",
                          "3h"  = "/{YYYY}{MM}{DD}{HH}00.kavrRAD.03.h5",
                          )


fc_file_path <- "/ec/res4/hpcperm/nhd/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/dini_eps/2025072112/"
fc_file_template <- "dini_tp_{YYYY}{MM}{DD}12+003_mbr000.grib2"

ob_file_format = "hdf5"
lead_time <- 12
fc_file_opts <- grib_opts(param_find = setNames(list(use_grib_shortName("tp")), parameter))



#my_param_defs <- modify_param_def(
#
#                               decription = "Total precipitation",
#                               #grib_opts(param_find = setNames(list(use_grib_shortName("tp")), parameter)),
#                               grib = new_grib_param(name="tp"),
#                               #grib = new_grib_param(
#                               #                      name = list(r = "tirf",
#                               #                                  g = "tgrp",
#                               #                                  s = "tsnowp"
#                               #                                  ),
#                               #                      ),
#                               #func = function(r, g, s) r + g + s,
#                               accum = extract_numeric(pcp_accum_period)
#                               )
#



#ob_file_opts <- hdf5_opts(data_path = "/pcp/data1/data", odim = FALSE, meta = TRUE)

#ob_file_opts   <- list(data_path="/dataset1/data1/data",
#                        odim="TRUE",
#                        meta=TRUE,
#                        invert_data=FALSE
#                        )
#

ob_file_opts   <- list(data_path="/dataset1/data1/data", #=NULL,
                        odim="TRUE",
                        meta=TRUE,
                        invert_data=F
                        )




#veri_time <- "202507121200"
##veri_time <- "202507120000"
init_time <- "202507210000" #this is the time when forecast was initialized
veri_time <- strftime(strptime(init_time, "%Y%m%d%H") + (lead_time[[1]] * 3600), "%Y%m%d%H00") #this is when the observation is avail 
print(veri_time)
fc_dttm  <- init_time

ob_file_name <- generate_filenames(file_path     = ob_file_path,
                                   file_date     = veri_time,
                                   file_template = ob_file_template)
fc_file_name <- generate_filenames(file_path     = fc_file_path,
                                   file_date     = init_time,
                                    lead_time = lead_time,
                                   file_template = fc_file_template)
# 2. Define the verification domain by reading the observation domain once.
#    verify_spatial will use this to regrid the forecast automatically.
#precip_ob <- read_grid(file_name=ob_file, parameter="RR", dttm=veri_time, file_format="netcdf", file_format_opts = ob_file_opts)

print(ob_file_name)
#precip_ob <- harpIO::read_grid(
#  ob_file_name,"Pcp",
#  hdf5_opts = hdf5_opts(data_path = "/pcp/data1/data", odim = FALSE, meta = TRUE)
#)


precip_ob <- read_grid(ob_file_name,
                     parameter        = parameter,
                     is_forecast      = TRUE, #Setting it to TRUE otherwise is not reading the data...
                     dttm             = veri_time,
                     lead_time        = NULL,
                     file_format      = ob_file_format,
                     file_opts = ob_file_opts)
                     #hdf5_opts = ob_file_opts)


precip_fc <- read_grid(fc_file_name,
                       parameter= "tp", #parameter,
                       #is_forecast=F,
                       #lead_time=lead_time,
                       dttm = veri_time, #fc_dttm,
                       file_format = "grib")
                       #file_opts = NULL) #fc_file_opts)
#
dmi_example <- read_grid(ob_file_name,
                         "pcp",
                         file_format = "hdf5",
                         file_format_opts = harpIO:::hdf5_opts(data_path=NULL,
                                                               odim=T,
                                                               meta=T,
                                                               invert_data=F))
ob_domain <- get_domain(dmi_example)

#plot_field(dmi_example)
print("Read grid from ob")
#ob_domain <- get_domain(
#  generate_filenames(
#    file_path     = "/hpcperm/kmek/obs/INCAPlus_1h/inca/",
#    file_template = "{YYYY}/{MM}/{DD}/INCAPlus_1h_RR_ANA_{YYYY}{MM}{DD}{HH}00.nc",
#    file_date     = "2025071212"
#  ),
#  file_format = "netcdf"
#)

dini_opts <- list(lead_time = lead_time,
                 veri_time = veri_time,
                 date_times = fc_dttm,
                 units = "kg/m2",
                 fc_model="dini_eps")

read_dini <- function(file_name,parameter="tp",lead_time,date_times,file_opts,...)
{

fc_data <- read_grid(file_name,
                       parameter= "tp", 
                       file_format = "grib")

veri_date <- strftime(strptime(date_times, "%Y%m%d%H") + (lead_time[[1]] * 3600), "%Y%m%d%H00") 

fc_model <- dini_opts["fc_model"]
fc_tbl <- tibble::tibble(
               valid_dttm   = dini_opts["veri_time"], #veri_time,
               parameter    = "Total Precipitation", #param$fullname,
               lead_time    = dini_opts["lead_time"], #lead_time,
               fcdate       = dini_opts["date_times"], #date_times,
               units        = dini_opts["units"], #units, #"kg/m2",
               gridded_data         =  geolist(fc_data)
               #!!as.name(fc_model) := geolist(fc_data)
)


return(fc_tbl)
}

print(fc_file_name)
test <- read_dini(file_name=fc_file_name,lead_time = lead_time, parameter="tp",date_times=fc_dttm,file_opts=dini_opts)

 # 3. Call verify_spatial with the correct arguments
slx_verification <- verify_spatial(
  dttm               =  fc_dttm, # veri_time, #"202507121200",
  parameter          = parameter, #"accrr1h",#"pcp",
  fcst_model         = "dini_eps",
  lead_time          = lead_time, #12 ,
  scores             = "SLX",
  window_sizes       = c(0, 5, 10, 15, 20,25,30), # c(1, 3, 5, 9, 17, 33), # These are the scales for SLX
  fc_file_path       = fc_file_path,
  fc_file_template   = fc_file_template,
  fc_file_format     = "dini", #"grib",
  fc_file_opts       = fc_file_opts,
  #fc_param_defs     = my_param_defs,
  fc_interp_method   = "nearest", # Method to regrid forecast to obs grid
  ob_interp_method   = NULL,
  ob_file_path       = ob_file_path,
  ob_file_template   = ob_file_template,
  ob_file_format     = ob_file_format,
  ob_file_opts       = ob_file_opts,
  ob_accumulation    = "3h", #"3h",
  fc_accumulation    = "3h",
  verif_domain       = ob_domain,
  return_data        = TRUE # To get the results back in a variable
)

# 4. Print the results
print(slx_verification$SLX)

