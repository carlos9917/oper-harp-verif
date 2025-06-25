

calcFSS2d <- function (members,observation,method,scales){
  # calculates the Member-Observations Agreement
  # members is list of two dimentionals geofield 
  # TODO: make sure that the number of dimensions is 2
   
  
  nmemb <- length(members)
  
  if (nmemb == 0 ){
    stop("calcFSS2d: zero members passed")
  }
  
  obnames <- names(observation)
  
  
  dims <- dim(observation[[1]])
  
  
  resdf <- NULL

  
  for (threperc in names(members[[1]])){
    

    summ <-  matrix(0.0, dims[1], dims[2])
    
    for( mIdx in seq(1,nmemb))
    {
      summ <- summ + members[[mIdx]][[threperc]]
    }
    
    for(scale in scales){

      fcprobability <- calcProbability(summ,scale)/nmemb
      
      obsprobability <- calcProbability(observation[[threperc]],scale)
      
      
      if (method == "MetOffice" ){
        fss <- fssVector(fcprobability,obsprobability)
      }
      else{
        if (method == "SMHI" ){
          fss <- fssSMHIVector(fcprobability,obsprobability)
        }
        else{
          stop(paste0("calcMOAgreement: method ",method, " is not supported."))
        }
      }
      
      if (is.null(resdf)){
        resdf <- tibble(label = threperc, scale=scale , fss = fss)
      }else{
        resdf <- add_row(resdf,label = threperc, scale=scale , fss = fss)
        
      }
      
    }
    
    
  
  }

  return(resdf)
  
}

