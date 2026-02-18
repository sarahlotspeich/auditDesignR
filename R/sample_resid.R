#' Sample residual sample audit based on naive model from Phase I variables.
#' @name sample_resid
#' @param formula Model formula used to calculate residuals, passed to \code{glm()} or \code{glm.nb()}.
#' @param family Type of model to be used to fit \code{formula}, passed to \code{glm()} or \code{glm.nb()}. Currently, the function accepts \code{family = "gaussian", "binomial", "poisson", "log-binomial",} or \code{"negbin"}.
#' @param dat Dataframe or matrix containing variables from \code{formula}.
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param wave1_Validated (For use with multi-wave designs) A logical vector with the validation indicator from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
#' @importFrom MASS glm.nb
#' @importFrom stats glm
sample_resid <- function(formula, family, dat, phI, phII, wave1_Validated = NULL) {
  ## If single wave, set wave1_validated = FALSE for all
  if (length(wave1_Validated) == 0) {
    wave1_Validated = rep(FALSE, phI)
  }
  
  ## Fit user-specified model
  if (family == "negbin") {
    fit = glm.nb(formula = as.formula(formula), 
                 data = dat)  
  } else if (family == "logbinom") {
    fit = glm(formula = as.formula(formula),
              family = binomial(link = "log"),
              data = dat)
  } else {
    fit = glm(formula = as.formula(formula), 
              family = family, 
              data = dat)
  }
  
  ## Get linear predictor
  mu = predict(fit)
  
  ## Calculate residuals (by model type)
  if(family == "gaussian") {
    ## Residuals linear regression: Y - mu
    r = as.vector(fit$y - mu)
  } else if (family == "binomial") {
    ## Residuals for logistic regression: Y - 1 / (1 + e ^ mu)
    r = as.vector(fit$y - (1 + exp(- mu)) ^ (- 1))
  } else if (family == "log-binomial") {
    ## Residuals for log-binomial regression: Y - N * (e ^ mu)
    y = fit$model[[1]][, 1] ### number of successes
    n = rowSums(fit$model[[1]]) ### number of trials
    r  = as.vector(y - n * exp(mu))
  } else if (family == "poisson") {
    ## Residuals for Poisson regression: Y - e ^ mu
    r = as.vector(fit$y - exp(mu))  
  } else if (family == "negbin") {
    ## Residuals for negative binomial regression (scaled): (Y - e ^ mu) / sqrt(e ^ mu - e ^ 2mu / theta)
    r = as.vector((fit$y - exp(mu)) / 
                    sqrt(exp(mu) - (exp(mu) ^ 2) / fit$theta))  
  }
  
  resid_id = data.frame(id = 1:phI, 
                        resid = r)
  
  ## If multi-wave, remove residuals from observations that were already validated
  resid_id = resid_id[!wave1_Validated, ]
  
  ## Order ascendingly by residuals
  resid_id = resid_id[order(resid_id$resid, decreasing = FALSE), ]
  ### Only validate smallest n/2 residuals
  smallest_resid = resid_id$id[1:(phII / 2)]
  
  ## Re-order descendingly by residuals
  resid_id = resid_id[order(resid_id$resid, decreasing = TRUE), ]
  ### Only validate smallest n/2 residuals
  largest_resid = resid_id$id[1:(phII / 2)]
  
  ## Return validation indicator
  V_resid = as.numeric(seq(1, phI) %in% c(smallest_resid, largest_resid))
  return(V_resid)
}
