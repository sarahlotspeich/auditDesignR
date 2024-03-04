#' Observed-data log-likelihood for measurement error settings with errors in outcome + exposure
#' @name od_loglik
#' @param params Vector of parameter vales for all models.
#' @param dat Data containing all columns referenced in \code{Y_unval}, \code{Y_val}, \code{X_unval}, \code{X_val}, \code{addl_covar}, and \code{Validated}.
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index). If outcome is error-free, \code{Y_unval = NULL} (DEFAULT).
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index). If covariates are error-free, \code{X_unval = NULL} (DEFAULT).
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param addl_covar Column(s) with additional error-free covariates (can be name or numeric index).
#' @param Validated Columns with the validation indicator (can be name or numeric index).
#' @param nondiff_Y_unval If \code{TRUE}, misclassification in \code{Y_unval} is assumed to be nondifferential, i.e., independent of \code{X_val} and \code{X_unval}. Default is \code{FALSE}.
#' @param nondiff_X_unval If \code{TRUE}, misclassification in \code{X_unval} is assumed to be nondifferential, i.e., independent of \code{Y_val}. Default is \code{FALSE}.
#' @param for_nlm If \code{TRUE} (DEFAULT), the function returns the negative function value (for use with built-in optimization routines like \code{nlm}).
#' @return Scalar function value.
#' @export
od_loglik <- function(params, dat, Y_val, Y_unval = NULL, X_val, X_unval = NULL, addl_covar = NULL, Validated, nondiff_Y_unval = FALSE, nondiff_X_unval = FALSE, for_nlm = TRUE) {
  N <- nrow(dat)
  dat[, "id"] <- seq(1, N)
  
  beta <- params[1:length(X_val)]
  eta <- params[-c(1:(length(X_val)))]
  
  # Parameters P(Y_val|X_val, addl_covar) - coeff for X_val were already saved as beta
  alpha <- eta[1:(1 + length(addl_covar))]
  eta <- eta[-c(1:(1 + length(addl_covar)))]
  
  # Parameters P(Y_unval|X_unval, Y_val, X_val, addl_covar)
  if (!is.null(Y_unval)) {
    if (!nondiff_Y_unval) {
      gamma_Ystar <- eta[1:(1 + length(c(Y_val, X_val, X_unval, addl_covar)))]
      eta <- eta[-c(1:(1 + length(c(Y_val, X_val, X_unval, addl_covar))))]
    } else {
      gamma_Ystar <- eta[1:(1 + length(c(Y_val, addl_covar)))]
      eta <- eta[-c(1:(1 + length(c(Y_val, addl_covar))))]
    }
  }
  
  # Parameters P(X_unval|Y_val, X_val, addl_covar)
  if (!is.null(X_unval)) {
    if (!nondiff_X_unval) {
      gamma_Xstar <- eta[1:(1 + length(c(Y_val, X_val, addl_covar)))]
      eta <- eta[-c(1:(1 + length(c(Y_val, X_val, addl_covar))))]
    } else {
      gamma_Xstar <- eta[1:(1 + length(c(X_val, addl_covar)))]
      eta <- eta[-c(1:(1 + length(c(X_val, addl_covar))))]
    }
  }
  
  # Parameters P(X_val|addl_covar)
  gamma_X <- eta[1:(1 + length(addl_covar))]
  
  # Validated subjects ----------------------------------------------------------------------
  val <- which(as.numeric(dat[, Validated]) == 1)
  n <- length(val)
  
  mu1 <- data.matrix(cbind(1, dat[val, c(addl_covar, X_val)])) %*% matrix(c(alpha, beta), ncol = 1)
  pY <- prob_logistic(y = dat[val, Y_val], mu = mu1)
  
  if (!is.null(Y_unval)) {
    if (!nondiff_Y_unval) {
      mu2 <- data.matrix(cbind(1, dat[val, c(X_unval, Y_val, X_val, addl_covar)])) %*% matrix(gamma_Ystar, ncol = 1)
    } else {
      mu2 <- data.matrix(cbind(1, dat[val, c(Y_val, addl_covar)])) %*% matrix(gamma_Ystar, ncol = 1)
    }
    pYstar <- prob_logistic(y = dat[val, Y_unval], mu = mu2)
  } else { pYstar <- rep(1, length(pY)) }
  
  if (!is.null(X_unval)) {
    if (!nondiff_X_unval) {
      mu3 <- data.matrix(cbind(1, dat[val, c(Y_val, X_val, addl_covar)])) %*% matrix(gamma_Xstar, ncol = 1)
    } else {
      mu3 <- data.matrix(cbind(1, dat[val, c(X_val, addl_covar)])) %*% matrix(gamma_Xstar, ncol = 1)
    }
    pXstar <- prob_logistic(y = dat[val, X_unval], mu = mu3)
  } else { pXstar <- rep(1, length(pY)) }
  
  mu4 <- data.matrix(cbind(1, dat[val, c(addl_covar)])) %*% matrix(gamma_X, ncol = 1)
  pX <- prob_logistic(y = dat[val, X_val], mu = mu4)
  
  l <- sum(log(pYstar) + log(pXstar) + log(pY) + log(pX))
  
  # Unvalidated subjects ---------------------------------------------------------------------
  unval <- which(as.numeric(dat[, Validated]) == 0)
  
  if (length(unval) > 0) {
    # Create complete dataset ------------------------------------------------------------------
    if (!is.null(X_unval)) {
      X_val_unique <- data.frame(unique(dat[val, X_val]))
      m <- nrow(X_val_unique)
      cd_unval <- dat[rep(unval, each = m), ]
      cd_unval[, X_val] <- X_val_unique[rep(seq(1, m), times = (N - n)), ]
    } else {
      cd_unval <- dat[unval, ]
      m <- 1
    }
    
    if (!is.null(Y_unval)) {
      cd_unval <- rbind(cd_unval, cd_unval)
      cd_unval[, Y_val] <- rep(c(0, 1), each = ((N - n) * m))
    }
    
    mu1 <- data.matrix(cbind(1, cd_unval[, c(addl_covar, X_val)])) %*% matrix(c(alpha, beta), ncol = 1)
    pY <- prob_logistic(y = cd_unval[, Y_val], mu = mu1)
    
    if (!is.null(Y_unval)) {
      if (!nondiff_Y_unval) {
        mu2 <- data.matrix(cbind(1, cd_unval[, c(X_unval, Y_val, X_val, addl_covar)])) %*% matrix(gamma_Ystar, ncol = 1)
      } else {
        mu2 <- data.matrix(cbind(1, cd_unval[, c(Y_val, addl_covar)])) %*% matrix(gamma_Ystar, ncol = 1)
      }
      pYstar <- prob_logistic(y = cd_unval[, Y_unval], mu = mu2)
    } else { pYstar <- rep(1, length(pY)) }
    
    if (!is.null(X_unval)) {
      if (!nondiff_X_unval) {
        mu3 <- data.matrix(cbind(1, cd_unval[, c(Y_val, X_val, addl_covar)])) %*% matrix(gamma_Xstar, ncol = 1)
      } else {
        mu3 <- data.matrix(cbind(1, cd_unval[, c(X_val, addl_covar)])) %*% matrix(gamma_Xstar, ncol = 1)
      }
      pXstar <- prob_logistic(y = cd_unval[, X_unval], mu = mu3)
    } else { pXstar <- rep(1, length(pY)) }
    
    mu4 <- data.matrix(cbind(1, cd_unval[, c(addl_covar)])) %*% matrix(gamma_X, ncol = 1)
    pX <- prob_logistic(y = cd_unval[, X_val], mu = mu4)
    
    l <- l + sum(log(rowsum(x = pY * pYstar * pXstar * pX, group = cd_unval[, "id"])))
  }
  
  if (for_nlm) { return(- l) } else { return(l) }
}
