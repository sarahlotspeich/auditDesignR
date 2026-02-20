#' Extreme-tail sampling based on a single Phase I variable.
#' @name sample_ets
#' @param ets_dat Vector containing one variable on which to sample. 
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param tails Which tails to sample from. Default is \code{tails = "both"}, but \code{"lower"} and \code{"upper"} only are possible.
#' @param wave1_Validated (For use with multi-wave designs) A logical vector with the validation indicator from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_ets <- function(ets_dat, phI, phII, tails = "both", wave1_Validated = NULL) {
  ## If single wave, set wave1_validated = FALSE for all
  if (length(wave1_Validated) == 0) {
    wave1_Validated = rep(FALSE, phI)
  }
  
  var_dat = data.frame(row_num = 1:phI, 
                       x = ets_dat) 
  
  if (tails == "both") {
    ## Order ascendingly by var_dat
    var_dat = var_dat[order(var_dat$x, decreasing = FALSE), ]
    
    ### Only validate smallest n/2 var_dat
    smallest_x = var_dat$row_num[1:(phII / 2)]
    
    ## Re-order descendingly by var_dat
    var_dat = var_dat[order(var_dat$x, decreasing = TRUE), ]
    
    ### Only validate smallest n/2 residuals
    largest_x = var_dat$row_num[1:(phII / 2)]
    
    ## Define validation indicator
    V_ets = as.numeric(1:phI %in% c(smallest_x, largest_x))
  } else if (tails == "lower") {
    ## Order ascendingly by var_dat
    var_dat = var_dat[order(var_dat$x, decreasing = FALSE), ]
    
    ### Only validate smallest n/2 var_dat
    smallest_x = var_dat$row_num[1:phII]
    
    ## Define validation indicator
    V_ets = as.numeric(1:phI %in% smallest_x)
  } else if (tails == "upper") {
    ## Re-order descendingly by var_dat
    var_dat = var_dat[order(var_dat$x, decreasing = TRUE), ]
    
    ### Only validate smallest n/2 residuals
    largest_x = var_dat$row_num[1:phII]
    
    ## Define validation indicator
    V_ets = as.numeric(1:phI %in% largest_x)
  }
  
  ## Return validation indicators
  return(V_ets)
}
