#' Prepares the observations to be an input for the spatial verification. It performs the Interpolation and converting the field into binary field if required.
#' @param start_date Date of the first forecast to read.
#' @param end_date Date of the last forecast to read.
#' @param by The time between forecasts. Should be a string of a number followed
#'   by a letter (the defualt is "6h"), where the letter gives the units - may
#'   be d for days, h for hours or m for minutes.
#' @param parameter The parameters to read as a character vector.
#' @param lead_time The lead times to read as a numeric vector.
#' @param lt_unit the unit of lead time, default hour "h"
#' @param ob_file_path The top level path for the forecast files to read.
#' @param ob_file_template The file type to generate the template for. Can be
#'   "harmoneps_grib", "harmeoneps_grib_fp", "harmoneps_grib_sfx", "meps_met",
#'   "harmonie_grib", "harmonie_grib_fp", "harmone_grib_sfx", "vfld", "vobs", or
#'   "fctable". If anything else is passed, it is returned unmodified. In this
#'   case substitutions can be used. Available substitutions are {YYYY} for
#'   year, \{MM\} for 2 digit month with leading zero, \{M\} for month with no
#'   leading zero, and similarly \{DD\} or \{D\} for day, \{HH\} or \{H\} for
#'   hour, \{mm\} or \{m\} for minute. Also \{LDTx\} for lead time and \{MBRx\}
#'   for ensemble member where x is the length of the string including leading
#'   zeros - can be omitted or 2, 3 or 4. Note that the full path to the file
#'   will always be file_path/template.
#' @param ob_file_format The format of the files to read. Can be e.g. "hdf5" or "grib".
#' @param ob_options A list with format-specific options for the reader function.
#' @param ob_interp_method Interpolation method to be used when transforming a forecast
#'   field to the verification grid.
#' @param ob_accumulation The accumulation type of the observation (or reference). This is only used for
#'   accumulated parameters (e.g. precipitation). NULL signifies that the field is accumulated
#'   from the start of the model run. That is probably rare for observations.
#'   Otherwise this should be a string containing a numerical value
#'   and a time unit, e.g. "15m" or "1h".
#' @param verif_domain A \code{geodomain} that defines the common verification grid.
#' @param nonbinary keep the parameter values as they are.
#' @param thresholds A vector of thresholds to convert the the parameter into binary fields.
#' @param percentiles A vector of percentiles to convert the parameter into binary fields.
#' @param save_to_dir Path where to save the processed values. If its value is NULL saving the output into files will be ignored.
#' @param returnData  I true then will return a list contains the processed values.
#' 
#' 
#' @return list of observations interpolated to the verif_domain. fields are accumulated based on discussion with Heiner (Nigel way)
#' 
#' @export

ens_spatial_calc_obs <- function(start_date,
                                 end_date,
                                 by                   = "6h",
                                 parameter,
                                 lead_time            = seq(0, 36, 6),
                                 lt_unit              = "h",
                                 ob_file_path         = NULL,
                                 ob_file_template     = NULL,
                                 ob_file_format       = "hdf5",
                                 ob_options           = list(),
                                 ob_interp_method     = "closest",
                                 ob_accumulation      = "1h",
                                 verif_domain         = NULL,
                                 nonbinary            = TRUE,
                                 thresholds           = NULL,
                                 percentiles          = NULL,                                 
                                 save_to_dir          = NULL,
                                 returnData           = TRUE
                                 
                                ) {
  
  
  if ( is.null(save_to_dir) & ( returnData == FALSE)){
    stop(" save_to_dir  is NULL and returnData is FALSE: No data to return. Will stop")
  }
  
  isThresholds = ! is.null(thresholds)
  isPercentiles = ! is.null(percentiles)
  

  if ( ! ( nonbinary | isThresholds |  isPercentiles )  ){
    stop("Nothing to do with the selected values of nonbinary, thresholds and percentiles. Will stop")
    
  }
  
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
  
  #message(paste(all_fc_dates, collapse = " "))
  #message(paste(all_ob_dates, collapse = " "))
  #message(paste(lead_time, collapse = " "))
  
  
  # if (return_data) {
  #   function_output <- list()
  #   list_counter    <- 0
  # }
  
  init <- list()
  
  get_ob <- function(obdate) {
    obfile <- get_filenames(
      file_date     = format(obdate, "%Y%m%d%H"),
      file_path     = ob_file_path,
      file_template = ob_file_template,
      parameter     = parameter
    )
    
    
    message("reading ", obfile)
    
    do.call(harpIO::read_grid,
            c(
              list(
                filename = obfile,
                file_format = ob_file_format,
                parameter = parameter
              ),
              ob_options
            ))
  }
  
  
  if (returnData){
    output <- list()
  }
  else{
    output <- NULL
  }
  
  ################# First crop the observations to the considered verif_domain   
  for (ob in seq_along(all_ob_dates)) {  # (obdate in all_ob_dates) looses POSIXct class
    
    obdate <- all_ob_dates[ob]
 
    strobdate <- format(obdate, "%Y%m%d%H%M") 

    message("=====\nobdate:", strobdate)
    
    obfield <- get_ob(obdate)
    if (inherits(obfield, "try-error")) { # e.g. missing observation
      if (harpenv$verbose) cat("Observation not found. Skipping.\n")
      next
    }
    if (prm$accum > 0) { # an accumulated field like precipitation
      if (is.null(ob_accumulation) || ob_accumulation < 0) {
        # RARE: observation is an accumulated field (e.g. reference run)
        # TODO: This does not look very useful, unless get_ob() can somehow deal with it.
        warning("Accumulated observation fields not yet validated.", immediate.=TRUE)
        obfield <- obfield - get_ob(obdate - prm$accum)
      } else {
        ostep <- readr::parse_number(ob_accumulation) * harpIO:::units_multiplier(ob_accumulation)
        if (ostep == prm$accum) { # this is easy !
          # nothing to do
        } else if  (ostep > prm$accum) { # this is easy !
          stop("The chosen accumulation time is smaller than that of the observations!")
        } else {
          
          nstep <- prm$accum / ostep
          for (i in 1:(nstep-1)) {
            obfield <- obfield + get_ob(obdate - i*ostep)
          }
        }
      }
    }

    # convert to common verification grid
    if (!is.null(ob_interp_method)) {
      if (is.null(init$regrid_ob)) {
        message("Initialising ob regridding.")
        init$regrid_ob <- meteogrid::regrid.init(
          olddomain = obfield,
          newdomain = verif_domain,
          method    = ob_interp_method)
      }
      obfield <- meteogrid::regrid(obfield, weights = init$regrid_ob)
      
      
    }
    
    ################# Do the Domain Accumulation Process  Without averaging.
    # Loop over my_matrix
    # fcfield[1, 1] remains the same
 

    data <- list()
    
    if (nonbinary){
      data[['nonbinary']] <- accumulateField(obfield)
    }
    
    
    attributes(data[['nonbinary']]) <- attributes(obfield)
    
    if (isPercentiles){
      
      percentileRes <-list()
      percentileThresholds = quantile(obfield,percentiles/100.)
      
      
      for (k in 1:length(percentiles)){
        precIndex <- format(percentiles[k], nsmall = 3)
        percentileRes[[precIndex]] <- accumulateBinField(obfield,percentileThresholds[k])
        attributes(percentileRes[[precIndex]]) <- attributes(obfield)
      }
      data[['percentile']] <- percentileRes
    }
    
     
    if (isThresholds){
      thresholdRes <-list()
      for (threshold in thresholds) {
        thrIndex <- format(threshold, nsmall = 3)
        thresholdRes[[thrIndex]] <- accumulateBinField(obfield,threshold)
        attributes(thresholdRes[[thrIndex]]) <- attributes(obfield)
      }
      data[['threshold']] <- thresholdRes
    }
    
    outfileBase <- paste0("grid_obs_", parameter, "_", as.character(strobdate))
    
    if ( !is.null(save_to_dir)){
       save_to_file <- file.path(save_to_dir,outfileBase)
       saveRDS(data, file = save_to_file)
    }
    
    if (returnData){
      output[[strobdate]] <- data
    }
    else{
      if (is.null(output))
      {
        output <- hashmap(strobdate, outfileBase)
      }
      else{
        output$insert(strobdate, outfileBase)
      }
    }
    
  }
 print(class(output))
 return(output)
}
