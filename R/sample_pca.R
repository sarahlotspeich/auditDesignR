#' Sample ETS-PCA audit based on Phase I variables.
#' @name sample_pca
#' @param pca_dat Dataframe or matrix containing variables on which to sample, passed to \code{princomp()}. This calculation is done using the correlation matrix, so these data do not need to be supplied on the scale. 
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param wave1_Validated (For use with multi-wave designs) A logical vector with the validation indicator from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_pca <- function(pca_dat, phI, phII, wave1_Validated = NULL) {
  ## If single wave, set wave1_validated = FALSE for all
  if (length(wave1_Validated) == 0) {
    wave1_Validated = rep(FALSE, phI)
  }
  
  ## Fit principal components analysis (PCA) on Phase I variables 
  pc = princomp(pca_dat, cor = TRUE)
  pc1 = data.frame(row_num = 1:phI, 
                   pc = pc$scores[, 1]) ### extract the first principal component
  
  ## If multi-wave, remove PCA from observations that were already validated
  pc1 = pc1[!wave1_Validated, ]
  
  ## Order ascendingly by first principal component
  pc1 = pc1[order(pc1$pc, decreasing = FALSE), ]
  ### Only validate smallest n/2 principal components
  smallest_pc1 = pc1$row_num[1:(phII / 2)]
  
  ## Re-order descendingly by first principal component
  pc1 = pc1[order(pc1$pc, decreasing = TRUE), ]
  ### Only validate smallest n/2 residuals
  largest_pc1 = pc1$row_num[1:(phII / 2)]

  ## Return validation indicator
  V_pca = as.numeric(seq(1, phI) %in% c(smallest_pc1, largest_pc1))
  return(V_pca)
}
