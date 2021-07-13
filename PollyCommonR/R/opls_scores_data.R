#' opls_scores_data
#' 
#' @param opl The "opls" object (ropls library) with typeC = 'OPLS-DA' for which scores are required.
#' 
#' @return A dataframe of scores and orthogonal component.
#' 
#' @examples 
#' data(scaurine)
#' attach(sacurine)
#' w = opls(dataMatrix,sampleMetadata[,'gender'],predI=1,orthoI=2)
#' opls_scores_data(w)
#' 
#' @import ropls
#' @export


opls_scores_data = function(opl){
  
  if(class(opl)!="opls"){
    stop("Object not of 'opls' class.")
  }
  if(opl@typeC!="OPLS-DA"){
    stop("You are not using a valid 'OPLS-DA' computed on the 'opls' object. Check the 'typeC' subclass in the opls s4 object.")
  }
  if(is.null(opl)){
    stop("Scores cannot be computed on NULL object.")
  }
  
  ortho = opl@orthoScoreMN #Model orthoScore
  scores = opl@scoreMN     #ModelScore
  comps = rbind(t(scores),t(ortho)) #binding both
  comps = as.data.frame(t(comps))  # as dataframe
  
  names(comps) = c("Score",paste("orthoScore",1:(ncol(comps)-1),sep=''))  # Naming columns in dataframe
  
  return(comps)
}