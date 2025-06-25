unlistpurrr <- function(inlist){
  
  
  k=1
  returnedValue <- list()
  for(i in 1:length(inlist)){
    subres <- inlist[[i]]
    for(j in 1:length(subres)){
      returnedValue[[k]] <- subres[[j]]
      k=k+1
    }
  }
  
  return(returnedValue)
}