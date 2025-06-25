#' @export 
ens_spatial_load_simulation <- function(ens_dates, 
                                        parameter            = "Accpcp3h",
                                        read_from_dir        = NULL,
                                        filesOnly=TRUE){
   
  
  init <- list()
  prm <- harpIO::parse_harp_parameter(parameter)

  output <- list()
  
  
  res <- purrr::pmap(ens_dates, function(lead_time,validdate,fcdate,nested){
    
    purrr::pmap(nested, function(eps_model,sub_model,actual_fcdate,member,lag_seconds,actual_lead_time_hours){
      
      result <- list()

      outfileBase <- paste0("grid_model_", parameter, "_", eps_model,"_", sub_model, "_mbr", formatC(member,width = 3, flag = "0"), "_", unixtime_to_str_datetime(fcdate,YMDhm), "+",formatC(lead_time,width = 4, flag = "0"))
      save_to_file <- file.path(read_from_dir, outfileBase)
      
      if (filesOnly == FALSE) {
        output[[outfileBase]] <- readRDS(result, file = save_to_file)
      }
      output[[outfileBase]]$fileBase <- outfileBase 
      
      
      return(output)
      
    })
    
  })
  
  
  returnedValue <- list()
  for(i in 1:length(res)){
    subres <- res[[i]]
    print(length(subres))
    for(j in 1:length(subres)){
      dat <- subres[[j]]
      
      for( k in names(dat)){
        print(k)
        returnedValue[[k]] <- dat[[k]]
      }
    }
  }
  return(returnedValue)
  
  
}

