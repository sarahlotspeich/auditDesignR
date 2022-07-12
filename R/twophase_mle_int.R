#' Maximum likelihood estimator for logistic regression with errors in binary outcome + exposure (with interaction)
#' @name twophase_mle_interact
#' @param dat Data containing all columns referenced in \code{Y_unval}, \code{Y_val}, \code{X_unval}, \code{X_val}, \code{addl_covar}, and \code{Validated}.
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index). If outcome is error-free, \code{Y_unval = NULL} (DEFAULT).
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index). If covariates are error-free, \code{X_unval = NULL} (DEFAULT).
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param addl_covar Column(s) with additional error-free covariates (can be name or numeric index). Default is \code{NULL}.
#' @param interact Column(s) with interaction term(s) (can be name or numeric index). Default is \code{NULL}.
#' @param Validated Columns with the validation indicator (can be name or numeric index).
#' @param int_Y_unval If \code{TRUE}, misclassification model for \code{Y_unval} includes \code{interact}. Default is \code{FALSE}.
#' @param int_X_unval If \code{TRUE}, misclassification model for \code{X_unval} includes \code{interact}. Default is \code{FALSE}.
#' @param noSE If \code{TRUE}, standard errors are calculated by inverting the hessian at convergence. Default is \code{FALSE}.
#' @return List of dataframes with model coefficients.
#' \item{mod_Y_unval}{Outcome misclassification model coefficients}
#' \item{mod_X_unval}{Exposure misclassification model coefficients}
#' \item{mod_Y_val}{Analysis model coefficients}
#' \item{mod_X_val}{Exposure model coefficients}
#' \item{conv}{\code{TRUE}/\code{FALSE} for convergence}
#' @export
twophase_mle_int <- function(dat, Y_val, Y_unval = NULL, X_val, X_unval = NULL, addl_covar = NULL, interact = NULL, Validated, int_Y_unval = FALSE, int_X_unval = FALSE, noSE = FALSE) {
  
  if (!is.null(Y_unval)) {
    if (int_Y_unval) {
      params_Y_unval <- c("(Intercept)", X_unval, Y_val, X_val, addl_covar, interact)
    } else {
      params_Y_unval <- c("(Intercept)", X_unval, Y_val, X_val, addl_covar) 
    }
  } else {
    params_Y_unval <- NULL
  }
  
  if (!is.null(X_unval)) {
    if (int_X_unval) {
      params_X_unval <- c("(Intercept)", Y_val, X_val, addl_covar, interact)  
    } else {
      params_X_unval <- c("(Intercept)", Y_val, X_val, addl_covar)
    }
  } else {
    params_X_unval <- NULL
  }
  
  params_Y_val <- c("(Intercept)", addl_covar, X_val)
  params_X_val <- c("(Intercept)", addl_covar)
  
  suppressWarnings(
    fit <- nlm(f = od_loglik_w_int,
               p = rep(0, length(c(params_Y_unval, params_X_unval, params_Y_val, params_X_val))),
               dat = dat,
               Y_val = Y_val,
               Y_unval = Y_unval,
               X_val = X_val,
               X_unval = X_unval,
               addl_covar = addl_covar,
               interact = interact,
               Validated = Validated,
               int_Y_unval = int_Y_unval, 
               int_X_unval = int_X_unval,
               iterlim = 250,
               steptol = 1E-4,
               Validated = Validated)
  )
  
  SE <- tryCatch(expr = sqrt(diag(solve(fit))),  error = function(err) {NA})
  conv <- fit$code <= 2
  if (conv) {
    beta <- fit$estimate[1]
    SE_beta <- SE[1]
    eta <- fit$estimate[-1]
    SE_eta <- SE[-1]
    alpha <- eta[1:(1 + length(addl_covar))]
    SE_alpha <- SE_eta[1:(1 + length(addl_covar))]
    eta <- eta[-c(1:(1 + length(addl_covar)))]
    SE_eta <- SE_eta[-c(1:(1 + length(addl_covar)))]
    
    coeff_Y_val <- data.frame(Coeff = params_Y_val,
                              Est = c(alpha[1], beta, alpha[-1]),
                              SE = c(SE_alpha[1], SE_beta, SE_alpha[-1]),
                              stringsAsFactors = FALSE)
    
    # Parameters P(Y_unval|X_unval, Y_val, X_val, addl_covar)
    if (!is.null(Y_unval)) {
      coeff_Y_unval <- data.frame(Coeff = params_Y_unval,
                                  Est = eta[1:length(params_Y_unval)],
                                  SE = SE_eta[1:length(params_Y_unval)],
                                  stringsAsFactors = FALSE)
      eta <- eta[-c(1:length(params_Y_unval))]
      SE_eta <- SE_eta[-c(1:length(params_Y_unval))]
    } else {
      coeff_Y_unval <- NULL
    }
    
    # Parameters P(X_unval|Y_val, X_val, addl_covar)
    if (!is.null(X_unval)) {
      coeff_X_unval <- data.frame(Coeff = params_X_unval,
                                  Est = eta[1:length(params_X_unval)],
                                  SE = SE_eta[1:length(params_X_unval)],
                                  stringsAsFactors = FALSE)
      eta <- eta[-c(1:length(params_X_unval))]
      SE_eta <- SE_eta[-c(1:length(params_X_unval))]
    } else {
      coeff_X_unval <- NULL
    }
    
    # Parameters P(X_val|addl_covar)
    coeff_X_val <- data.frame(Coeff = params_X_val,
                              Est = eta[1:length(params_X_val)],
                              SE = SE_eta[1:length(params_X_val)],
                              stringsAsFactors = FALSE)
    
    return(list(mod_Y_val = coeff_Y_val, mod_Y_unval = coeff_Y_unval, mod_X_unval = coeff_X_unval, mod_X_val = coeff_X_val, conv = conv))
    
  } else {
    return(list(mod_Y_val = NA, mod_Y_unval = NA, mod_X_unval = NA, mod_X_val = NA, conv = conv))
  }
}
