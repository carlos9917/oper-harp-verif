#' Prepares the model output to be an input for the spatial verification. It performs the Interpolation and converting the field into binary field if required.
#' @param ens_dates output of ens_spatial_prep_dates
#' @param end_date Date of the last forecast to read.
#' @param parameter The parameters to read as a character vector.
#' @param fc_file_path The top level path for the forecast files to read.
#' @param fc_file_template The file type to generate the template for. Can be
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
#' @param fc_file_format The format of the files to read. Can be e.g. "fa" or "grib".
#' @param fc_options A list with format-specific options for the reader function.
#' @param fc_interp_method Interpolation method to be used when transforming a forecast
#'   field to the verification grid.
#' @param fc_accumulation The accumulation type of the forecast. This is only used for
#'   accumulated parameters (e.g. precipitation). NULL signifies that the field is accumulated
#'   from the start of the model run. Otherwise this should be a string containing a numerical value
#'   and a time unit, e.g. "15m" or "1h".
#' @param verif_domain A \code{geodomain} that defines the common verification grid.
#' @param nonbinary keep the parameter values as they are.
#' @param thresholds A vector of thresholds to convert the the parameter into binary fields.
#' @param percentiles A vector of percentiles to convert the parameter into binary fields.
#' @param save_to_dir Path where to save the processed values. If its value is NULL saving the output into files will be ignored.
#' @param returnData  I true then will return a list contains the processed values.
#' 
#' 
#' @return list of values interpolated to the verif_domain. fields are accumulated based on discussion with Heiner (Nigel way)
#' @export
ens_spatial_calc_simulation <- function(ens_dates, 
                                        parameter            = "Accpcp3h",
                                        fc_file_path         = NULL,
                                        fc_file_template     = NULL,
                                        fc_file_format       = "fa",
                                        fc_options           = list(),
                                        fc_interp_method     = "closest",
                                        fc_accumulation      = NULL,                         
                                        verif_domain         = NULL,
                                        nonbinary            = TRUE,
                                        thresholds           = NULL,
                                        percentiles          = NULL,                                        
                                        save_to_dir          = NULL,
                                        returnData           = TRUE
){
  
  if ( is.null(save_to_dir) & ( returnData == FALSE)){
    stop(" save_to_dir  is NULL and returnData is FALSE: No data to return. Will stop")
  }
  
  isThresholds = ! is.null(thresholds)
  isPercentiles = ! is.null(percentiles)
  

  if ( ! ( nonbinary | isThresholds |  isPercentiles )  ){
    stop("Nothing to do with the selected values of nonbinary, thresholds and percentiles.")
    
  }
  
  
  init <- list()
  
  
  prm <- harpIO::parse_harp_parameter(parameter)
  # TODO: Check the accumulation units well
  prm$accum <- round(prm$accum /3600)
  
 
  print(prm)
  
  
  
  read_file_grid <- function(filePath) {
    message("reading ", filePath)
    do.call(harpIO::read_grid,
            c(list(filename = filePath, file_format = fc_file_format, 
                   parameter = parameter),
              fc_options))
  }
  
  
  get_filePath <- function(actual_fcdate, actual_lead_time_hours,eps_model,sub_model,member) {
    get_filenames(
      file_path  = fc_file_path,
      start_date = unixtime_to_str_datetime(actual_fcdate,YMDhm),
      end_date   = unixtime_to_str_datetime(actual_fcdate,YMDhm),
      lead_time  = actual_lead_time_hours,
      parameter  = parameter,
      eps_model  = eps_model,
      sub_model  = sub_model,
      members    = member,
      file_template = fc_file_template,
      filenames_only = TRUE)
    
  }
  
  
  if (returnData){
    output <- list()
  }
  
  res <- purrr::pmap(ens_dates, function(lead_time,validdate,fcdate,nested){
    
    purrr::pmap(nested, function(eps_model,sub_model,actual_fcdate,member,lag_seconds,actual_lead_time_hours){
      file_path = get_filePath(actual_fcdate    =  unix2datetime(actual_fcdate) ,
                               actual_lead_time_hours = actual_lead_time_hours,
                               eps_model        = eps_model,
                               sub_model        = sub_model, 
                               member           = member)
      
      
      
      
      print(file_path)
      if(!file.exists(file_path)){
        message("file does not exist, Ignoring Date", unix2datetime(actual_fcdate), ", lead_time ", actual_lead_time_hours,", File Path", file_path)
      }
      else{
        
        fcfield <- read_file_grid(file_path)
        
        if (prm$accum > 0) {
          if (is.null(fc_accumulation) || fc_accumulation < 0) {
            #TODO: check actual_lead_time_hours units every where : beter to be inseconds
            if (actual_lead_time_hours > prm$accum) { # if ldt==accum, you don't need to decumulate
              file_path = get_filePath(actual_fcdate    =  unix2datetime(actual_fcdate) ,
                                       actual_lead_time_hours = actual_lead_time_hours - prm$accum,
                                       eps_model        = eps_model,
                                       sub_model        = sub_model,
                                       member           = member)
              if(file.exists(file_path)){
                fcfield <- fcfield - read_file_grid(file_path)
              }
              else
              {
                stop(pate0("ens_spatial_calc_simulation: file ", file_path, " does not exists"))
              }
            }
          } 
          else {
            # In rare cases the forecast model needs "accumulating" rather than "decumulating"
            #       e.g. when verifying INCA against radar
            fstep <- readr::parse_number(fc_accumulation) * 
              harpIO:::units_multiplier(fc_accumulation)
            if (fstep == prm$accum) { # this is easy !
              # nothing to do
            } else if  (fstep > prm$accum) { # this is easy !
              stop("The chosen accumulation time is smaller than what is available in the forecasts!")
            } else {
              nstep <- prm$accum / fstep
              for (i in 1:(nstep - 1)) {
                file_path = get_filePath(actual_fcdate    =  unix2datetime(actual_fcdate) ,
                                         actual_lead_time_hours = actual_lead_time_hours - i * fstep,
                                         eps_model        = eps_model,
                                         sub_model        = sub_model,
                                         member           = member)
                if(!file.exists(file_path)){
                  fcfield <- fcfield - read_file_grid(file_path)
                }
                else
                {
                  fcfield[,] = NA
                }
              }
            }
          }
        }
        
        
        # convert forecast to common verification grid
        if (!is.null(fc_interp_method)) {
          if (is.null(init$regrid_fc)) {
            message("Initialising fc regridding.")
            init$regrid_fc <- 
              meteogrid::regrid.init(
                olddomain = fcfield,
                newdomain = verif_domain,
                method = fc_interp_method
              )
          }
          fcfield <- meteogrid::regrid(fcfield, weights = init$regrid_fc)
          
          ################# Do the Domain Accumulation Process  Without averaging.
          # Loop over my_matrix
          # fcfield[1, 1] remains the same
          
          result <- list() 
          
          if (nonbinary){
            result[['nonbinary']] <- accumulateField(fcfield)
            attributes(result[['nonbinary']]) <- attributes(fcfield)
          }
          
          
          if (isPercentiles){
            
            percentileRes <-list()
            percentileThresholds = quantile(fcfield,percentiles/100.)
            
            for (k in 1:length(percentiles)) {
              precIndex <- format(percentiles[k], nsmall = 3)
              percentileRes[[precIndex]] <- accumulateBinField(fcfield,percentileThresholds[k])
              attributes(percentileRes[[precIndex]]) <- attributes(fcfield)
              
              
            }
            
            result[['percentile']] <- percentileRes
            
          }
          
          
          if (isThresholds){
            thresholdRes <-list()
            for (threshold in thresholds) {
              thrIndex <- format(threshold, nsmall = 3)
              thresholdRes[[thrIndex]] <- accumulateBinField(fcfield,threshold)
              attributes(thresholdRes[[thrIndex]]) <- attributes(fcfield)
            }
            result[['threshold']] <-thresholdRes
          }
          
          ####################### End OF Domain Accumulation Process################
          
          #lead_time,validdate,fcdate,eps_model,sub_model,member
          outfileBase <- paste0("grid_model_", parameter, "_", eps_model,"_", sub_model, "_mbr", formatC(member,width = 3, flag = "0"), "_", unixtime_to_str_datetime(fcdate,YMDhm), "+",formatC(lead_time,width = 4, flag = "0"))
          result$fileBase <- outfileBase
          if (! is.null(save_to_dir)){
            save_to_file <- file.path(save_to_dir, outfileBase)
            saveRDS(result, file = save_to_file)
          }
          if (returnData){
            output[[outfileBase]] <- result
          }
          
        }
      }
      if (returnData){
        return(output)
      }      
    })
    
  })
  
  if (returnData){
    returnedValue <- list()
    for(i in 1:length(res)){
      subres <- res[[i]]
      for(j in 1:length(subres)){
        dat <- subres[[j]]
        
        for( k in names(dat)){
          returnedValue[[k]] <- dat[[k]]
        }
      }
    }
    return(returnedValue)
  }
  
  
}

