#' Observed-data log-likelihood for measurement error settings with errors in outcome + exposure
#' @name od_loglik
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index).
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index).
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param Validated Columns with the validation indicator (can be name or numeric index).
#' @return Scalar function value.
#' @export
od_loglik <- function(params, dat, Y_val, Y_unval, X_val, X_unval, Validated) {
  beta <- params[1]
  eta <- params[-1]

  alpha <- eta[1]
  gamma_Ystar <- eta[2:5]
  gamma_Xstar <- eta[6:8]
  gamma_X <- eta[9]

  val <- which(as.numeric(dat[, Validated]) == 1)

  # Validated subjects ----------------------------------------------------------------------
  mu1 <- data.matrix(cbind(1, dat[val, X_val])) %*% matrix(c(alpha, beta), ncol = 1)
  pY <- prob_logistic(y = dat[val, Y_val], mu = mu1)

  mu2 <- data.matrix(cbind(1, dat[val, c(X_unval, Y_val, X_val)])) %*% matrix(gamma_Ystar, ncol = 1)
  pYstar <- prob_logistic(y = dat[val, Y_unval], mu = mu2)

  mu3 <- data.matrix(cbind(1, dat[val, c(Y_val, X_val)])) %*% matrix(gamma_Xstar, ncol = 1)
  pXstar <- prob_logistic(y = dat[val, X_unval], mu = mu3)

  mu4 <- matrix(rep(gamma_X, length(val)), ncol = 1)
  pX <- prob_logistic(y = dat[val, X_val], mu = mu4)

  l <- sum(log(pYstar) + log(pXstar) + log(pY) + log(pX))

  # Unvalidated subjects ---------------------------------------------------------------------
  cd_unval <- rbind(dat[-val, ], dat[-val, ])
  cd_unval[, Y_val] <- rep(c(0, 1), each = nrow(dat[-val, ]))
  cd_unval <- rbind(cd_unval, cd_unval)
  cd_unval[, X_val] <- rep(c(0, 1), each = 2 * nrow(dat[-val, ]))

  mu1 <- data.matrix(cbind(1, cd_unval[, X_val])) %*% matrix(c(alpha, beta), ncol = 1)
  pY <- prob_logistic(y = cd_unval[, Y_val], mu = mu1)

  mu2 <- data.matrix(cbind(1, cd_unval[, c(X_unval, Y_val, X_val)])) %*% matrix(gamma_Ystar, ncol = 1)
  pYstar <- prob_logistic(y = cd_unval[, Y_unval], mu = mu2)

  mu3 <- data.matrix(cbind(1, cd_unval[, c(Y_val, X_val)])) %*% matrix(gamma_Xstar, ncol = 1)
  pXstar <- prob_logistic(y = cd_unval[, X_unval], mu = mu3)

  mu4 <- matrix(rep(gamma_X, nrow(cd_unval)), ncol = 1)
  pX <- prob_logistic(y = cd_unval[, X_val], mu = mu4)

  l <- l + sum(log(rowsum(x = pY * pYstar * pXstar * pX, group = rep(seq(1, nrow(dat[-val, ])), times = 4))))

  return(-l)
}
