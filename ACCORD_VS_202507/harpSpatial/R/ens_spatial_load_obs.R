#' @export
 
ens_spatial_load_obs <- function(start_date,
                                 end_date,
                                 by                   = "6h",
                                 parameter,
                                 lead_time            = seq(0, 36, 6),
                                 lt_unit              = "h",
                                 read_from_dir          = NULL
                                ) {

  
  # TODO: we may need more options! masked interpolation, options by score,
  prm <- harpIO::parse_harp_parameter(parameter)
  by_secs <- harpIO:::units_multiplier(by) * readr::parse_number(by)
  
  
  
  # Loop over dates to prevent excessive data volumes in memory
  
  pctStart_date <-
    lubridate::ymd_hm(start_date, tz = "UTC", truncated = 4)
  pctEnd_date <-
    lubridate::ymd_hm(end_date, tz = "UTC", truncated = 4)
  
  
  all_fc_dates <-
    seq(as.numeric(pctStart_date) ,
        as.numeric(pctEnd_date),
        by_secs) %>%  lubridate::as_datetime()
  
  
  lt_scale <- harpIO:::units_multiplier(lt_unit)
  
  #lead_time <- lead_time * lt_scale on Alex Version
  
  all_ob_dates <-
    (rep(all_fc_dates, each = length(lead_time)) + lead_time * lt_scale)   %>%
    unique() %>%
    sort()
  
  

    output <- list()

  
  ################# First crop the observations to the considered verif_domain   
  for (ob in seq_along(all_ob_dates)) {  # (obdate in all_ob_dates) looses POSIXct class
    
    obdate <- all_ob_dates[ob]
 
    strobdate <- format(obdate, "%Y%m%d%H%M") 
    
    outfileBase <- paste0("grid_obs_", parameter, "_", as.character(strobdate))
    
    output[[strobdate]] <- readRDS(file.path(read_from_dir,outfileBase))
    output[[strobdate]]$fileBase <- outfileBase

  }
  

 return(output)
}
