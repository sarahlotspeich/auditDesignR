#' Sample Fisher's information sample audit based on naive model from Phase I variables.
#' @name sample_fisher
#' @param formula Model formula used to calculate residuals, passed to \code{glm()} or \code{glm.nb()}.
#' @param family Type of model to be used to fit \code{formula}, passed to \code{glm()} or \code{glm.nb()}. Currently, the function accepts \code{family = "gaussian", "binomial", "poisson", "log-binomial",} or \code{"negbin"}.
#' @param dat Dataframe or matrix containing variables from \code{formula}.
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param tails Which tails to sample from. Default is \code{tails = "both"}, but \code{"lower"} and \code{"upper"} only are possible.
#' @param wave1_Validated (For use with multi-wave designs) A logical vector with the validation indicator from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
#' @importFrom MASS glm.nb
#' @importFrom stats glm
sample_fisher <- function(formula, family, dat, phI, phII, tails = "both", wave1_Validated = NULL) {
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
  
  ## Get model matrix 
  mod_mat = model.matrix(fit)
  
  ## Calculate deviances (by model type)
  if(family == "gaussian") {
    ## Mean of X
    xbar = mod_mat[, 2]
    ## Fisher information contribution for linear regression: (Y - mu) ^ 2 / sigma^2
    fi = as.vector((mod_mat[, 2] - xbar) ^ 2 / sigma(fit) ^ 2)
  } else if (family == "binomial") {
    ## Predicted probability of event
    p = (1 + exp(- mu)) ^ (- 1)
    ## Fisher information weight
    w = p * (1 - p)
    ## Weighted mean of X = sum(w * X) / sum(w)
    xbar_w = sum(w * mod_mat[, 2]) / sum(w)
    ## Fisher Information contributions for logistic regression
    fi = as.vector(w * (mod_mat[, 2] - xbar_w) ^ 2) 
  } else if (family == "poisson") {
    ## Predicted mean count
    lambda = exp(mu)
    ## Calculate mean predicted count = sum(mu * x) / sum(mu)
    xbar_mu = sum(mu * mod_mat[, 2]) / sum(mu)
    ## Fisher Information contributions for Poisson regression
    fi = as.vector(mu * (mod_mat[, 2] - xbar_mu) ^ 2) 
  } else if (family == "log-binomial") {
    ## Predicted mean count
    lambda = exp(mu)
    ## Fisher information weight
    w = lambda / (1 - lambda)
    ## Calculate mean predicted count = sum(w * x) / sum(w)
    xbar_w = sum(w * mod_mat[, 2]) / sum(w)
    ## Fisher Information contributions for log-binomial regression
    fi = as.vector(w * (mod_mat[, 2] - xbar_w) ^ 2) 
  } else if (family == "negbin") {
    ## Predict means 
    lambda = exp(mu) 
    ## Fisher Information contributions for negative binomial regression
    fi = as.vector(lambda * mod_mat[, 2] ^ 2 / (1 + fit$theta * lambda)) 
  }
  
  fi_id = data.frame(id = 1:phI, 
                     fi = fi)
  
  ## If multi-wave, remove Fisher informations from observations that were already validated
  fi_id = fi_id[!wave1_Validated, ]
  
  if (tails == "both") {
    ## Order ascendingly by deviances
    fi_id = fi_id[order(fi_id$fi, decreasing = FALSE), ]
    ### Only validate smallest n/2 devuals
    smallest_fi = fi_id$id[1:(phII / 2)]
    
    ## Re-order descendingly by deviances
    fi_id = fi_id[order(fi_id$fi, decreasing = TRUE), ]
    ### Only validate smallest n/2 deviances
    largest_fi = fi_id$id[1:(phII / 2)]
    
    ## Return validation indicator
    V_fi = as.numeric(seq(1, phI) %in% c(smallest_fi, largest_fi))
  } else if (tails == "lower") {
    ## Order ascendingly by deviances
    fi_id = fi_id[order(fi_id$fi, decreasing = FALSE), ]
    ### Only validate smallest n/2 devuals
    smallest_fi = fi_id$id[1:phII]
    ## Return validation indicator
    V_fi = as.numeric(seq(1, phI) %in% smallest_fi)
  } else if (tails == "upper") {
    ## Re-order descendingly by deviances
    fi_id = fi_id[order(fi_id$fi, decreasing = TRUE), ]
    ### Only validate smallest n/2 deviances
    largest_fi = fi_id$id[1:phII]
    ## Return validation indicator
    V_fi = as.numeric(seq(1, phI) %in% largest_fi)
  }
  return(V_fi)
}
