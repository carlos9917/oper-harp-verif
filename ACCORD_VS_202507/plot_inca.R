library(harp)
library(rlang)
library(here)

#filename <- "/hpcperm/kmek/obs/INCAPlus_1h/inca/2025/07/01/INCAPlus_1h_RR_ANA_202507010400.nc"
veri_time <- "202507011600"
veri_time <- "202507121200"
parameter <- "acc1h"
parameter <- "accrr1h"

ob_file_path <- "/hpcperm/kmek/obs/INCAPlus_1h/inca"

ob_file_format = "netcdf"
ob_file_template <- "{YYYY}/{MM}/{DD}/INCAPlus_1h_RR_ANA_{YYYY}{MM}{DD}{HH}00.nc"

# for do nicer plotting:
observation <- "incaPlus"
ob_units <- "mm/h"

plot_path <- "./PLOTS/"
ob_plot_name <- "IncaPlus_test.png"

script_path <- "."
plotting_functions_file <- paste0(script_path, "/plotting_functions.R")


ob_dttm         <- veri_time
ob_dttm_POSIXct <- as.POSIXct(veri_time, "%Y%m%d%H", tz = "UTC")

### read observation ###
ob_file_name <- generate_filenames(file_path     = ob_file_path,
                                   file_date     = ob_dttm,
                                   file_template = ob_file_template)
ob_file_opts     <- netcdf_opts(proj4_var  = "lambert_conformal_conic",
                             #x_rev      = TRUE,
                             #y_rev      = TRUE,
                             param_find = list2(!!parameter := "RR") # TODO: make this flexible depending on the accumulation period.
)

ob_data <- read_grid(ob_file_name,
                     parameter        = parameter,
                     is_forecast      = FALSE,
                     dttm             = ob_dttm,
                     file_format      = ob_file_format,
                     file_format_opts = ob_file_opts)



### quick harp plotting ###
#plot_field(ob_data)

### advanced plotting ###

# create a tibble
ob_tbl <- tibble::tibble(
               valid_dttm   = ob_dttm,
               parameter    = parameter,
               # lead_time    = lead_time,
               fcdate       = ob_dttm,
               units        = ob_units,
               !!as.name(observation) := geolist(ob_data)
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

ob_title <- paste(parameter, "incaPlus", ob_dttm_POSIXct)

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

