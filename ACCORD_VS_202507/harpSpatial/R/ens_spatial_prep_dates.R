#' Prepares the dates for each ensemble ensemble member.
#' @param start_date Date of the first forecast to be read in. Should be in
#'   YYYYMMDDhh format. Can be numeric or charcter.
#' @param end_date Date of the last forecast to be read in. Should be in
#'   YYYYMMDDhh format. Can be numeric or charcter.
#' @param by The time between forecasts. Should be a string of a number followed
#'   by a letter (the defualt is "6h"), where the letter gives the units - may
#'   be d for days, h for hours or m for minutes.
#' @param lead_time The lead times to read as a numeric vector.
#' @param lt_unit the unit of lead time, default hour "h"
#' @param ens_members the output of ens_spatial_prep_members
#' 
#' @return A nested tibble with colums leadtime, validdate, fcstdate.
#' 
#' @export 
ens_spatial_prep_dates <- function(start_date,
                                   end_date,
                                   by                   = "6h",
                                   lead_time            = seq(0, 36, 6),
                                   lt_unit              = "h",
                                   ens_members          = NULL
){
  
  # ens_members <- tidyr::unnest(ens_members, tidyr::one_of("eps_model")) 
  # ens_members <- tidyr::unnest(ens_members, tidyr::one_of("sub_model")) 
  # ens_members <- tidyr::unnest(ens_members, tidyr::one_of("member","members_out","lag"))
  #  
  # 
  # for(i in 1:nrow(ens_members)) {
  #   row <- ens_members[i,]
  #   print(row)
  # }
  # 
  # res <- ens_members %>%  purrrlyr::by_row(print)
  # return(0)
  
 
 
  # TODO: we may need more options! masked interpolation, options by score,
  by_secs <- harpIO:::units_multiplier(by) * readr::parse_number(by)

  all_dates <- seq_dates(start_date, end_date, by)
  
  
  lt_scale <- harpIO:::units_multiplier(lt_unit)
  
  units_multiplier_vec <- Vectorize(harpIO:::units_multiplier, USE.NAMES = FALSE)
  
  output <- NULL
  
  
  ##TODO: if (tidyr_new_interface()) {
  if(TRUE){
    ens_unnested_members <- tidyr::unnest(ens_members, -tidyr::one_of("eps_model"))
  } else {
    ens_unnested_members <- tidyr::unnest(ens_members)
  }
  
  
  for (fcst_date in all_dates) {
    
    data_files <- ens_unnested_members %>%
      dplyr::transmute(
        file_names = purrr::pmap(
          list(
            eps_model = .data$eps_model,
            sub_model = .data$sub_model,
            members   = .data$member,
            lags      = .data$lag
          ),
          function(eps_model, sub_model, members, lags) get_ens_dates(
            start_date     = fcst_date,
            end_date       = fcst_date,
            by             = by,
            lags           = lags,
            eps_model      = eps_model,
            sub_model      = sub_model,
            lead_time      = lead_time,
            members        = members
          )
          
          
          
        )
      )
    
 
  
  
    
    ##TODO:L if (tidyr_new_interface()) {
    if(TRUE){
      data_files  <- tidyr::unnest(data_files, tidyr::one_of("file_names"))
      #ens_members <- tidyr::unnest(ens_members, -tidyr::one_of("eps_model"))
    } else {
      data_files  <- tidyr::unnest(data_files)
      #ens_members <- tidyr::unnest(ens_members)
    }


 # data_files <- dplyr::left_join(
 #  data_files,
 #  ens_members,
 #  by = c("eps_model", "sub_model", "member")
 # )
 # #   Get the data
    
    message("Preparing date: ", fcst_date)
    
    #read_function <- get(paste("read", file_format, "interpolate", sep = "_"))
    forecast_data <- data_files %>%
      dplyr::mutate(
        ####Str actual_fcdateStr = unixtime_to_str_datetime(.data$actual_fcdate,, YMDhm),
        ####Str member = paste0("mbr", formatC(.data$member, width = 3, flag = "0")),        
        member = .data$member,
        validdate = .data$actual_fcdate + .data$lead_time * 3600 + .data$lag_seconds, 
        ####Str validdateStr = unixtime_to_str_datetime(validdate, YMDhm),
        fcdate = actual_fcdate + .data$lag_seconds,
        actual_lead_time_hours = .data$lead_time + .data$lag_seconds/3600 # TODO: make sure that the lag and lead time are in hours
        
      )         
                    

   if (is.null(output)){
      output <- forecast_data
    }
    else{
     output <- dplyr::bind_rows(output,forecast_data)    
    }
    
 }

output <- output %>% tidyr::nest(nested = - c(lead_time,validdate,fcdate))
  
return(  output )
  


}