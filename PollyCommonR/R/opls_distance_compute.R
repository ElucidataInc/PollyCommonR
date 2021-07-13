#' opls_distance_compute
#' 
#' @param opl The "opls" object (ropls library) with typeC = 'OPLS-DA'.
#' 
#' @param significance A boolean indicating whether significant points are to be computed.
#' @param level A float indicating the significance level to be used in the two-sided chi-square test. Ignored is signiicance=FALSE.
#' 
#' @return If significance=TRUE, A list of 1) dataframe containg scores distance data, 
#' 2) a vector containing distances corresponding to significance level, and, 
#' 3) A vector of indices of significant points in the dataframe.
#' Otherwise only (1) is returned.
#' 
#' @import stats
#' @export
#' 
#' 


opls_distance_compute =  function(opl, significance = TRUE, level = 0.05){
  
  #### code copied and modified from plot modality of ropls library ####
  
  tCompMN <- opl@scoreMN
  pCompMN <- opl@loadingMN
  
  
  tCompMN <- cbind(opl@scoreMN[, 1], opl@orthoScoreMN[, 1])
  pCompMN <- cbind(opl@loadingMN[, 1], opl@orthoLoadingMN[, 1])
  colnames(pCompMN) <- colnames(tCompMN) <- c("h1", paste("o", 1, sep = ""))
  
  mahInvCovMN = solve(cov(tCompMN))
  
  pcaResMN <- cbind(sdsVn = apply(tCompMN,
                                  1,
                                  function(x) sqrt(t(as.matrix(x)) %*% mahInvCovMN %*% as.matrix(x))),
                    dsVn = apply(opl@suppLs[["xModelMN"]] - tcrossprod(tCompMN, pCompMN),
                                 1,
                                 function(x) sqrt(drop(crossprod(x[complete.cases(x)])))))
  
  
  ########################## copied till here ###########################    
  
  pcaResMN = as.data.frame(pcaResMN)   
  
  if (significance==TRUE){
    if(class(level)=='numeric'){
      if(level <1 || level>0){
        pcaResThrVn <- c(sqrt(qchisq(1-level/2, 2)),
                         (mean(pcaResMN[, 2]^(2/3)) + sd(pcaResMN[, 2]^(2/3)) * qnorm(0.975))^(3/2))
        
        pcaResExtVi <- union(which(pcaResMN[, 1] > pcaResThrVn[1]),
                             which(pcaResMN[, 2] > pcaResThrVn[2]))
        
        return (list(pcaResMN, pcaResThrVn, pcaResExtVi))
        
      }
      
      else{
        warning("Use a valid significance level. Ignoring significance argument.")
      }
    }
    else{
      warning("Significance level has to be 'numeric'. Ignoring significance argument.")
    }
  }
  return(pcaResMN)
  
}


