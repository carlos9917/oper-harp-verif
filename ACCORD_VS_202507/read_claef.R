library(harp)
library(rlang)
parameter <- "acc1h"

filename <- "/hpcperm/kmek/models/CLAEF1k/20250707/00/CLAEF00+0006:00.grb2"
filename <- "/hpcperm/kmek/models/CLAEF1k/20250701/00/CLAEF02+0016:00.grb"
filename <- "/hpcperm/kmek/models/CLAEF1k/20250701/00/CLAEF00+0002:00.grb"

parameter <- "pcp"
fc_file_opts     <- grib_opts( param_find = setNames( list(list(key = 'indicatorOfParameter', value = 61)),
                              parameter))  # not needed if parameter="pcp" - needed to make it uniform with INCA 

#fc_file_opts     <- grib_opts( param_find = setNames( list(list(key = 'indicatorOfParameter', value = 61),list(key = 'typeOfLevel', value = "surface")),
#                              parameter))  # not needed if parameter="pcp" - needed to make it uniform with INCA 

#fc <- read_grid(file_name=filename, parameter=parameter, dttm="20250701", file_format="grib", file_format_opts = fc_file_opts)


fc_file_opts     <- grib_opts( param_find = setNames( list(list(key = 'indicatorOfParameter', value = 61)), "pcp"))
precip_fc <- read_grid(file_name=filename, parameter=parameter, dttm="2025070102", file_format="grib", file_format_opts = fc_file_opts)
#plot_field(fc)
