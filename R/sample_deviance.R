#' Sample deviance sample audit based on naive model from Phase I variables.
#' @name sample_deviance
#' @param formula Model formula used to calculate deviances, passed to \code{glm()}.
#' @param family Type of model to be used to fit \code{formula}, passed to \code{glm()}. Currently, the function accepts \code{family = "gaussian", "binomial",} or \code{"poisson"}.
#' @param dat Dataframe or matrix containing variables from \code{formula}.
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param wave1_Validated (For use with multi-wave designs) A logical vector with the validation indicator from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_deviance <- function(formula, family, dat, phI, phII, wave1_Validated = NULL) {
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
  
  ## Calculate deviances (by model type)
  if(family == "gaussian") {
    ## Deviance for linear regression: (Y - mu) ^ 2
    d = as.vector((fit$y - mu) ^ 2)
  } else if (family == "binomial") {
    ## Predicted probability of event 
    p = (1 + exp(- mu)) ^ (- 1)
    ## Deviance for logistic regression: 2 * Y - 1 / (1 + e ^ mu)
    d = 2 * (fit$y * log(fit$y / p) + (1 - fit$y) * log((1 - fit$y) / (p)))
  } else if (family == "poisson") {
    ## Predicted mean count
    lambda = exp(mu)
    ## Deviances for Poisson regression: Y log(Y / lambda) - (Y - lambda)
    d = as.vector(2 * (fit$y * log(fit$y / lambda) - (fit$y - lambda)))
  }
  
  dev_id = data.frame(id = 1:phI, 
                      dev = d)
  
  ## If multi-wave, remove deviances from observations that were already validated
  dev_id = dev_id[!wave1_Validated, ]
  
  ## Order ascendingly by deviances
  dev_id = dev_id[order(dev_id$dev, decreasing = FALSE), ]
  ### Only validate smallest n/2 devuals
  smallest_dev = dev_id$id[1:(phII / 2)]
  
  ## Re-order descendingly by deviances
  dev_id = dev_id[order(dev_id$dev, decreasing = TRUE), ]
  ### Only validate smallest n/2 deviances
  largest_dev = dev_id$id[1:(phII / 2)]
  
  ## Return validation indicator
  V_dev = as.numeric(seq(1, phI) %in% c(smallest_dev, largest_dev))
  return(V_dev)
}
