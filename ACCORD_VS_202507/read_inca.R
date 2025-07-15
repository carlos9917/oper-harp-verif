library(harp)
library(rlang)
parameter <- "acc1h"
filename <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/INCAPlus_1h/inca/2025/07/01/INCAPlus_1h_RR_ANA_202507010400.nc"
ob_file_opts     <- netcdf_opts(proj4_var  = "lambert_conformal_conic", param_find = list2(!!parameter := "RR"))
read_grid(file_name=filename, parameter="RR", dttm="202507010400", file_format="netcdf", file_format_opts = ob_file_opts)
