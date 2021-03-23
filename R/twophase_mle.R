#' Maxium likelihood estimator for logistic regression with errors in binary outcome + exposure
#' @name twophase_mle
#' @param dat Data containing all columns referenced in \code{Y_unval}, \code{Y_val}, \code{X_unval}, \code{X_val}, \code{addl_covar}, and \code{Validated}.
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index). If outcome is error-free, \code{Y_unval = NULL} (DEFAULT).
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index). If covariates are error-free, \code{X_unval = NULL} (DEFAULT).
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param addl_covar Column(s) with additional error-free covariates (can be name or numeric index).
#' @param Validated Columns with the validation indicator (can be name or numeric index).
#' @param nondiff_Y_unval If \code{TRUE}, misclassification in \code{Y_unval} is assumed to be nondifferential, i.e., independent of \code{X_val} and \code{X_unval}. Default is \code{FALSE}.
#' @param nondiff_X_unval If \code{TRUE}, misclassification in \code{X_unval} is assumed to be nondifferential, i.e., independent of \code{Y_val}. Default is \code{FALSE}.
#' @param noSE If \code{TRUE}, standard errors are calculated by inverting the hessian at convergence. Default is \code{FALSE}.
#' @return List of dataframes with model coefficients.
#' \item{mod_Y_unval}{Outcome misclassification model coefficients}
#' \item{mod_X_unval}{Exposure misclassification model coefficients}
#' \item{mod_Y_val}{Analysis model coefficients}
#' \item{mod_X_val}{Exposure model coefficients}
#' \item{conv}{\code{TRUE}/\code{FALSE} for convergence}
#' @export
twophase_mle <- function(dat, Y_val, Y_unval = NULL, X_val, X_unval = NULL, addl_covar = NULL, Validated, nondiff_Y_unval = FALSE, nondiff_X_unval = FALSE, noSE = FALSE) {
  if (!nondiff_Y_unval) {
    params_Y_unval <- c("(Intercept)", X_unval, Y_val, X_val, addl_covar)
  } else {
    params_Y_unval <- c("(Intercept)", Y_val, addl_covar)
  }

  if (!nondiff_X_unval) {
    params_X_unval <- c("(Intercept)", Y_val, X_val, addl_covar)
  } else {
    params_X_unval <- c("(Intercept)", X_val, addl_covar)
  }

  params_Y_val <- c("(Intercept)", addl_covar, X_val)
  params_X_val <- c("(Intercept)", addl_covar)

  suppressWarnings(
    fit <- nlm(f = od_loglik,
               p = rep(0, length(c(params_Y_unval, params_X_val, params_Y_val, params_X_val))),
               dat = dat,
               Y_val = Y_val,
               Y_unval = Y_unval,
               X_val = X_val,
               X_unval = X_unval,
               addl_covar = addl_covar,
               Validated = Validated,
               nondiff_X_unval = nondiff_X_unval,
               nondiff_Y_unval = nondiff_Y_unval,
               iterlim = 250,
               hessian = F)
  )

  conv <- fit$code <= 2
  if (conv) {

  } else {

  }
}
