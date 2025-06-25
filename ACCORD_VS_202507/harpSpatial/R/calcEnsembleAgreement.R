
calcMMAgreement <- function (members, method = "MetOffice", slim=5,alpha=0.5,epsilon=0.0001){
  # calculates the Member-Member Agreement
  # members is list of two dimentionals geofield 
  # TODO: make sure that the number of dimensions is 2
  
  nmemb <- length(members)
  combs <- combn(1:nmemb,2)
  cdims <- dim(combs)
  mdims <- dim(members[[1]])
  
  result <-  matrix(0.0, mdims[1], mdims[2])
  
  for( cbIdx in seq(1,cdims[2]))
  {
    if (method == "MetOffice" ){
      calcPrimaryAgreementMetOfficeV2(members[[combs[1,cbIdx]]],members[[combs[2,cbIdx]]],result,slim,alpha,epsilon)
    }
    else{
      if (method == "SMHI" ){
        calcPrimaryAgreementSMHIV2(members[[combs[1,cbIdx]]],members[[combs[2,cbIdx]]],result,slim,alpha,epsilon)
      }
      else{
        stop(paste0("calcMMAgreement: method ",method, " is not supported."))
      }
    }
  }
  
  
  result <- result/cdims[2]
  
  attributes(result) <- attributes(members[[1]])
  
  return(result)
}
 

calcMOAgreement <- function (members, observation,method = "MetOffice",  slim=5,alpha=0.5,epsilon=0.0001){
  # calculates the Member-Observations Agreement
  # members is list of two dimentionals geofield 
  # TODO: make sure that the number of dimensions is 2
  
  nmemb <- length(members)
  
  dims <- dim(observation)
  
  result <-  matrix(0.0, dims[1], dims[2])
  
  for( mIdx in seq(1,nmemb))
  {
    if (method == "MetOffice" ){
      calcPrimaryAgreementMetOfficeV2(members[[mIdx]],observation,result,slim,alpha,epsilon)
    }
    else{
      if (method == "SMHI" ){
        calcPrimaryAgreementSMHIV2(members[[mIdx]],observation,result,slim,alpha,epsilon)
      }
      else{
        stop(paste0("calcMOAgreement: method ",method, " is not supported."))
      }
    }
    
  }
  
  
  result <- result/nmemb
  
  attributes(result) <- attributes(members[[1]])
  
  return(result)
  
  
}



