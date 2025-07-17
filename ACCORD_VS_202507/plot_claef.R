library(here)
library(harp)
library(yaml)
library(rlang)

########################
#veri_time <- "202507010300"
#parameter <- "accrr3h"

veri_time <- "202507011600"
veri_time <- "202507121200"
parameter <- "acc1h"
parameter <- "pcp" #this works too for reading but breaks the plot_field
parameter <- "accrr1h"


fc_file_path <- "/hpcperm/kmek/models/CLAEF1k"

fc_file_format = "grib"
fc_file_template <- "{YYYY}{MM}{DD}/00/CLAEF00+0016:00.grb" #for the old files 
fc_file_template <- "{YYYY}{MM}{DD}/00/CLAEF01+0014:00.grb2" #for the new files
fc_file_template <- "{YYYY}{MM}{DD}/00/{det_model}+{LDT4}:00.grb2" #for the new files

#fc_file_path <- "/ment_arch3/aladin/PRECIP_ARCH"
#fc_file_template <- "{YYYY}{MM}{DD}/{det_model}+{LDT4}.grb"
lead_time <- 12

fc_file_format = "grib"

fcst_model <- "claef_1k_00"
fcst_model <- "CLAEF00"
# for do nicer plotting:
fc_units <- "mm/h"

plot_path <- "./PLOTS"
fc_plot_name <- "claef1k_grb_test.png"

script_path <- "."
plotting_functions_file <- paste0(script_path, "/plotting_functions.R")

########################
ob_dttm         <- veri_time
fc_dttm         <- strftime(strptime(veri_time, "%Y%m%d%H") - (lead_time * 3600), "%Y%m%d%H")
fc_dttm_POSIXct <- as.POSIXct(fc_dttm, "%Y%m%d%H", tz = "UTC")

param               <- parse_harp_parameter(parameter)

### read observation ###
fc_file_name <- generate_filenames(file_path     = fc_file_path,
                                   file_date     = fc_dttm,
                                   lead_time     = lead_time,
                                   file_template = fc_file_template,
                                   det_model     = fcst_model)

fc_file_opts     <- grib_opts(
                          param_find = setNames(
                                            list(list(key = 'indicatorOfParameter',
                                                      value = 61)),
			  parameter))  # not needed if parameter="pcp" - needed to make it uniform with INCA ##TODO: make flexible for different accumulation times


#this definition for the new grib2 files
fc_file_opts     <- grib_opts(
                          param_find = setNames(list(use_grib_shortName('tp')), parameter))



#                      dttm             = fc_dttm,
#                      fcst_model       = fcst_model,
#                      parameter        = parameter,
#                      lead_time        = lead_time,
#                      file_path        = fc_file_path,
#                      file_template    = fc_file_template,
#                      file_format      = fc_file_format,
#                      file_format_opts = fc_file_opts,
#                      param_defs       = fc_param_defs,
#                      return_data      = TRUE
# )

fc_data <- read_forecast(
                     dttm             = fc_dttm,
                     fcst_model       = fcst_model,
                     parameter        = parameter,
                     lead_time        = lead_time,
                     file_path        = fc_file_path,
                     file_template    = fc_file_template,
                     file_format      = fc_file_format,
                     file_format_opts = fc_file_opts,
                     param_defs       = list(),
                     return_data      = TRUE
)


fc_data <- read_grid(fc_file_name,
                     parameter        = parameter,
                     is_forecast      = TRUE,
                     dttm             = fc_dttm,
                     file_format      = fc_file_format,
                     file_format_opts = fc_file_opts,
		     param_defs       = list(),
                     lead_time        = lead_time)



### quick harp plotting ###
plot_field(fc_data)

### advanced plotting ###

# create a tibble
fc_tbl <- tibble::tibble(
               valid_dttm   = ob_dttm,
               parameter    = param$fullname,
               lead_time    = lead_time,
               fcdate       = fc_dttm,
               units        = fc_units,
               !!as.name(fcst_model) := geolist(fc_data)
)

# source functions for plotting #
source(plotting_functions_file)

# set breaks, define palette and title #
breaks    = c(0, 0., 1, 3, 5, 7, 10., 15.,
                 20, 25, 30, 40, 50, 60, 70, 80, 100)

palette   = c("#e5ebec", "#c7e5fb","#8cc7f2", "#45a6eb", "#1c73b2",
                 "#4d991b", "#71ce9c", "#b4db72", "#f5f305",
                 "#f6d125", "#f6a625", "#f54125", "#ae092f",
                 "#d59de5", "#9c04c6", "#23052b")
print(parameter)
print(fcst_model)
print(fc_dttm_POSIXct)
print(lead_time)
print(param$acc_unit)

fc_title <- paste0(parameter, "   ",
                fcst_model, "   ",
                fc_dttm_POSIXct, " + ",
                lead_time, param$acc_unit
)


if (!dir.exists(plot_path)){
	dir.create(plot_path, recursive = TRUE)
}

# plot field #
fc_gg <- plot_panel_field(fc_tbl,
                          fcst_model,
                          title     = fc_title,
                          breaks    = breaks,
                          palette   = palette,
                          plot_path = plot_path,
                          plot_name = fc_plot_name
)


