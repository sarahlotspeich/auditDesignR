#' Sample score function sample audit based on naive model from Phase I variables.
#' @name sample_sfs
#' @param formula Model formula used to calculate residuals, passed to \code{glm()}.
#' @param family Type of model to be used to fit \code{formula}, passed to \code{glm()}. Currently, the function accepts \code{family = "gaussian", "binomial",} or \code{"poisson"}.
#' @param dat Dataframe or matrix containing variables from \code{formula}.
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param X Column with the error-prone covariate incorporated into sampling (can be name or numeric index).
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_sfs <- function(formula, family, dat, phI, phII, X) {
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
  smallest_sf = which(order(sf, decreasing = FALSE) <= phII / 2) ## only validate smallest n/2 residuals
  largest_sf = which(order(sf, decreasing = TRUE) <= phII / 2) ## ... and largest n/2 

  V_sfs = as.numeric(seq(1, phI) %in% c(smallest_sf, largest_sf))
  return(V_sfs)
}
