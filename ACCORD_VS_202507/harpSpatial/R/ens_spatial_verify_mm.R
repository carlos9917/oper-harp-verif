#' Based on Dey et al. 2016
#' @export
ens_spatial_verify_mm <- function(ens_dates              = NULL,
                                  parameter              = NULL,
                                  read_simulation_from   = NULL, # may be directory or tible of  data framee passed
                                  slim                   = 5,
                                  alpha                  = 0.5,
                                  epsilon                = 0.01,
                                  save_to_dir            = NULL,
                                  returnData             = TRUE,
                                  isSMHI                 = FALSE
){
  
  
  

  if ( is.null(read_simulation_from)){
    stop("ens_spatial_verify_mm: read_simulation_from is NULL. Nothing to do. will stop")
  }
  
  if ( is.null(save_to_dir) & (!returnData)){
    stop("ens_spatial_verify_mm: save_to_dir is NULL and returnData is False. No place specified to save or return data. Will stop")
  }
  
  
  if ( is.character(read_simulation_from)){
    isLoadSimulation= TRUE
  }
  else{
    isLoadSimulation = FALSE
  }
  
  
  
  res = purrr::pmap(ens_dates, function(lead_time,validdate,fcdate,nested){
    
    output <- list()
    
    
    
    model_files <- purrr::pmap(nested, function(eps_model,sub_model,actual_fcdate,member,lag_seconds,actual_lead_time_hours){
      
      model_file <- paste0("grid_model_", parameter, "_", eps_model,"_", sub_model, "_mbr", formatC(member,width = 3, flag = "0"), "_", unixtime_to_str_datetime(fcdate,YMDhm), "+",formatC(lead_time,width = 4, flag = "0"))
      return(model_file) 
      
      
    })
    
    
    model_files <- unlist(model_files)
    
    fcdataNonbinary  <- list()
    
    for( model_file in model_files){
      if(isLoadSimulation){
        modelData <- readRDS(file.path(read_simulation_from,model_file))
      }
      else{
        modelData <-read_simulation_from[[model_file]]
      }
      fcdataNonbinary[[length(fcdataNonbinary) + 1]] <- modelData[["nonbinary"]]
    }
    
    validdateStr <- unixtime_to_str_datetime(validdate,YMDhm)
    
    
    
    fcdateStr <- unixtime_to_str_datetime(fcdate,YMDhm)
    lead_timeStr <- formatC(lead_time,width = 4, flag = "0")
    

    if (isSMHI){
      res <- calcMMAgreement(members      = fcdataNonbinary,
                                 method       = "SMHI",
                                 slim         = slim,
                                 alpha        = alpha,
                                 epsilon      = epsilon)
      
      outputfileBase <- paste0("mmsv_", parameter, "_",fcdateStr, "+",lead_timeStr)
    }
    else
    {
      res <- calcMMAgreement(members      = fcdataNonbinary,
                             method       = "MetOffice",
                             slim         = slim,
                             alpha        = alpha,
                             epsilon      = epsilon)
      
      outputfileBase <- paste0("mmsv_", parameter, "_",fcdateStr, "+",lead_timeStr)
    }

    
    
    
    if(!is.null(save_to_dir)){
      res_file <- file.path(save_to_dir,outputfileBase)
      saveRDS(res,res_file)
    }
    if(returnData){
      outputItem = list( file = res_file, validdate = validdate, fcdate = fcdate, lead_time = lead_time,fileName = outputfileBase, spatialVerificationResult = res)
      
    }
    else{
      outputItem = list( file = res_file, validdate = validdate, fcdate = fcdate, lead_time = lead_time,fileName = outputfileBase)
    }
    output[[length(output) + 1]] <- outputItem
    
    
    return(output)
    
  })
  
  
  k=1
  returnedValue <- list()
  for(i in 1:length(res)){
    subres <- res[[i]]
    for(j in 1:length(subres)){
      returnedValue[[k]] <- subres[[j]]
      k=k+1
    }
  }
  return(returnedValue)
  
  
  
}
