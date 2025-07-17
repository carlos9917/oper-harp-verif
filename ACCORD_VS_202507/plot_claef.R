library(harp)
library(rlang)
library(here)

#filename <- "/hpcperm/kmek/obs/INCAPlus_1h/inca/2025/07/01/INCAPlus_1h_RR_ANA_202507010400.nc"
veri_time <- "202507011600"
veri_time <- "202507121200"
parameter <- "acc1h"
parameter <- "accrr1h"
parameter <- "pcp"

fc_file_path <- "/hpcperm/kmek/models/CLAEF1k"

fc_file_format = "grib"
fc_file_template <- "{YYYY}{MM}{DD}/00/CLAEF00+0016:00.grb" #for the old files 
fc_file_template <- "{YYYY}{MM}{DD}/00/CLAEF01+0014:00.grb2" #for the new files
fc_file_template <- "{YYYY}{MM}{DD}/00/CLAEF00+0012:00.grb2" #for the new files

# for do nicer plotting:
observation <- "CLAEF1k"
fc_units <- "mm/h"

plot_path <- "./PLOTS/"
ob_plot_name <- "CLAEF_test.png"

script_path <- "."
plotting_functions_file <- paste0(script_path, "/plotting_functions.R")


fc_dttm         <- veri_time
fc_dttm_POSIXct <- as.POSIXct(veri_time, "%Y%m%d%H", tz = "UTC")

### read observation ###
fc_file_name <- generate_filenames(file_path     = fc_file_path,
                                   file_date     = fc_dttm,
                                   file_template = fc_file_template)
#This definition okay for old grib1 files
#fc_file_opts     <- grib_opts( param_find = setNames( list(list(key = 'indicatorOfParameter', value = 61)),
#                              parameter))  # not needed if parameter="pcp" - needed to make it uniform with INCA 


#this definition for the new grib2 files
fc_file_opts     <- grib_opts(
                          param_find = setNames(list(use_grib_shortName('tp')), parameter))


#fc <- read_grid(file_name=filename, parameter=parameter, dttm="2025070102", file_format="grib", file_format_opts = fc_file_opts)
#plot_field(fc)

fc_data <- read_grid(fc_file_name,
                     parameter        = parameter,
                     is_forecast      = TRUE, #FALSE,
                     dttm             = fc_dttm,
                     file_format      = fc_file_format,
                     lead_time = 12,
                     file_format_opts = fc_file_opts)

print(fc_data)

### quick harp plotting ###
#plot_field(fc_data)

### advanced plotting ###

# create a tibble
ob_tbl <- tibble::tibble(
               valid_dttm   = fc_dttm,
               parameter    = parameter,
               # lead_time    = lead_time,
               fcdate       = fc_dttm,
               units        = fc_units,
               !!as.name(observation) := geolist(fc_data)
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

ob_title <- paste(parameter, "CLAEF1k", fc_dttm_POSIXct)

if (!dir.exists(plot_path)){
        dir.create(plot_path, recursive = TRUE)
}
# plot field #
ob_gg <- plot_panel_field(ob_tbl,
                          observation,
                          title     = ob_title,
                          breaks    = breaks,
                          palette   = palette,
                          plot_path = plot_path,
                          plot_name = ob_plot_name
)

#ob_file_opts     <- netcdf_opts(proj4_var  = "lambert_conformal_conic", param_find = list2(!!parameter := "RR"))
#df <- read_grid(file_name=filename, parameter="RR", dttm=veri_time, file_format="netcdf", file_format_opts = ob_file_opts) #,data_frame=TRUE)
#
#plot_field(df)

