#' compute_opls_s_data
#' 
#' Computes the loading and correlation dataframe required for S-plot.
#' 
#' @param opl An opls object
#' @return A dataframe
#' @examples
#' compute_s_data(opls_object)
#' @import ropls stats
#' @export
compute_opls_s_data = function(opl=NULL){
  
  if(is.null(opl)){
    warning("NULL object cannot be used.")
    return(NULL)
  }
  if(class(opl)!="opls"){
    stop("Object not of 'opls' class.")
  }
  if(opl@typeC!="OPLS-DA"){
    warning("You are not using a valid OPLS-DA computed on the 'opls' object. Check the 'typeC' subclass in the opls s4 object.")
  }
  
  
  require("ropls")
  Data = as.data.frame(cbind(opl@loadingMN,opl@loadingMN/(sd(opl@scoreMN)*opl@xSdVn)))
  names(Data) = c("Loading","Correlation")
  
  return(Data)
  
}