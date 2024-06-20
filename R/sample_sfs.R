#' Sample score function sample audit based on naive model from Phase I variables.
#' @name sample_sfs
#' @param formula Model formula used to calculate residuals, passed to \code{glm()}.
#' @param family Type of model to be used to fit \code{formula}, passed to \code{glm()}. Currently, the function accepts \code{family = "gaussian", "binomial",} or \code{"poisson"}.
#' @param dat Dataframe or matrix containing variables from \code{formula}.
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param X Column with the error-prone covariate incorporated into sampling (can be name or numeric index).
#' @param absValue Logical, if \code{absValue = FALSE} (the default) the design is based on the weighted score functions, rather than their absolute values.
#' @param wave1_Validated (For use with multi-wave designs) A logical vector with the validation indicator from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_sfs <- function(formula, family, dat, phI, phII, X, absValue = FALSE, wave1_Validated = NULL) {
  ## If single wave, set wave1_validated = FALSE for all
  if (length(wave1_Validated) == 0) {
    wave1_Validated = rep(FALSE, phI)
  }
  
  ## Fit user-specified model
  fit = glm(formula = as.formula(formula), 
            family = family, 
            data = dat)
  
  ## Get linear predictor
  mu = predict(fit)
  
  ## Calculate residuals (by model type)
  if(family == "gaussian") {
    ## Residuals linear regression: Y - mu
    r = as.vector(fit$y - mu)
  } else if (family == "binomial") {
    ## Residuals for logistic regression: Y - 1 / (1 + e ^ mu)
    r = as.vector(fit$y - (1 + exp(- mu)) ^ (- 1))
  } else if (family == "poisson") {
    ## Residuals for Poisson regression: Y - e ^ mu
    r = as.vector(fit$y - exp(mu))  
  }
  
  ## Multiply residuals by X 
  sf = r * dat[, X]
  
  sf_id = data.frame(id = 1:phI, 
                     sf = sf)
  
  ## If multi-wave, remove residuals from observations that were already validated
  sf_id = sf_id[!wave1_Validated, ]
  
  ## Order ascendingly by residuals
  if (absValue) {
    sf_id = sf_id[order(abs(sf_id$sf), decreasing = FALSE), ]
  } else {
    sf_id = sf_id[order(sf_id$sf, decreasing = FALSE), ]
  }
  ### Only validate smallest n/2 residuals
  smallest_sf = sf_id$id[1:(phII / 2)]
  
  ## Re-order descendingly by residuals
  if(absValue) {
    sf_id = sf_id[order(abs(sf_id$sf), decreasing = TRUE), ] 
  } else {
    sf_id = sf_id[order(sf_id$sf, decreasing = TRUE), ]  
  }
  ### Only validate smallest n/2 residuals
  largest_sf = sf_id$id[1:(phII / 2)]
  
  ## Return validation indicator
  V_sfs = as.numeric(seq(1, phI) %in% c(smallest_sf, largest_sf))
  return(V_sfs)
}
