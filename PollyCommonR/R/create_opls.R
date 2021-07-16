#' create_opls
#' 
#' Creates a class opls object for OPLS-DA
#' 
#' @param sample_data Input dataframe with shape (n x m) where n is the number of samples
#' and m = number of variables
#' @param metadata Sample metadata dataframe with number of rows = n
#' @param condition A column in metadata on which the OPLS-DA has to be performed. Should
#' have binary classes.
#' @param scalC Character: either no centering nor scaling ('none'), mean-centering only ('center'),
#' mean-centering and pareto scaling ('pareto'), or mean-centering and unit variance scaling ('standard') (default)
#' @returns An 'opls' class object (see library 'ropls')
#' @examples  
#' create_opls(sample_data, metadata, condiition)
#' @import ropls
#' @export
create_opls = function(sample_data,metadata,condition,
                       scalC = c("none", "center", "pareto", "standard")[4]){
  
  message("Create OPLS object started...")
  
  require(ropls)
  
  if (nrow(sample_data) == 0){
    warning("Not a valid input matrix, have NANs or infs in the matrix")
    return (NULL)
  }
  
  if (nrow(metadata) == 0){
    warning("Input metadata is not a valid matrix, have NANs or infs in the matrix")
    return (NULL)
  }
  
  if(nrow(sample_data)!=nrow(metadata)){
    stop("Input sample_data and metadata must have the same number of rows.")
  }
  
  if(!(condition %in% names(metadata))){
    stop("Input condition is not a valid column in metadata.")
  }
  
  if(!(length(unique(metadata[,condition]))==2)){
    stop("OPLS-DA is only available for binary targets. Select a different condition and try again.")
  }
  
  if(!(scalC %in% c("none", "center", "pareto", "standard"))){
    warning("Invalid scaling option selected. Using default: 'standard'.")
    scalc = "standard"
  }
  
  message("Making the opls object...")
  crossval = min(nrow(sample_data),7)
  
  opls_obj = opls(sample_data,metadata[,condition], scaleC = scalC, predI = 1, orthoI = 1, fig.pdfC = FALSE)
  
  message("opls object returned.")
  return(opls_obj)
  
}