#' Score vectors for measurement error setting with errors in outcome + exposure
#' @name score
#' @param comp_dat Dataframe with complete data.
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index).
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index).
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param Validated Columns with the validation indicator (can be name or numeric index).
#' @param beta Parameter value for the log odds ratio of \code{X_val} on \code{Y_val}.
#' @param eta Parameter values for other nuisance parameters.
#' @return Matrix with one column for the score vector for each parameter \code{beta}, \code{eta}.
#' @export
score <- function(comp_dat, Y_val, Y_unval, X_val, X_unval, Validated, beta, eta) {
  val <- which(comp_dat[, Validated] == 1)
  n <- length(val)

  alpha <- eta[1]
  gamma_Ystar <- eta[2:5]
  gamma_Xstar <- eta[6:8]
  gamma_X <- eta[9]

  Si_theta <- data.frame(
    comp_dat,
    Si_beta = NA, Si_alpha = NA,
    Si_gamma0 = NA, Si_gamma1 = NA, Si_gamma2 = NA,
    Si_gamma3 = NA, Si_gamma4 = NA, Si_gamma5 = NA,
    Si_gamma6 = NA, Si_gamma7 = NA, row.names = NULL)

  ## S(beta), S(alpha)
  mu1 <- data.matrix(cbind(1, Si_theta[, X_val])) %*% matrix(c(alpha, beta), ncol = 1)
  pY <- prob_logistic(y = Si_theta[, Y_val], mu = mu1)
  d_alpha <- (2 * Si_theta[, Y_val] - 1) * sigmoid(mu1) * (1 - sigmoid(mu1))
  d_beta <- Si_theta[, X_val] * d_alpha

  ## S(gamma0), ... , S(gamma3)
  mu2 <- data.matrix(cbind(1, Si_theta[, c(X_unval, Y_val, X_val)])) %*% matrix(gamma_Ystar, ncol = 1)
  pYstar <- prob_logistic(y = Si_theta[, Y_unval], mu = mu2)
  d_gamma0 <- (2 * Si_theta[, Y_unval] - 1) * sigmoid(mu2) * (1 - sigmoid(mu2))
  d_gamma1 <- Si_theta[, X_unval] * d_gamma0
  d_gamma2 <- Si_theta[, Y_val] * d_gamma0
  d_gamma3 <- Si_theta[, X_val] * d_gamma0

  ## S(gamma4), ..., S(gamma6)
  mu3 <- data.matrix(cbind(1, Si_theta[, c(Y_val, X_val)])) %*% matrix(gamma_Xstar, ncol = 1)
  pXstar <- prob_logistic(y = Si_theta[, X_unval], mu = mu3)
  d_gamma4 <- (2 * Si_theta[, X_unval] - 1) * sigmoid(mu3) * (1 - sigmoid(mu3))
  d_gamma5 <- Si_theta[, Y_val] * d_gamma4
  d_gamma6 <- Si_theta[, X_val] * d_gamma4

  ## S(gamma7)
  mu4 <- matrix(rep(gamma_X, nrow(Si_theta)), ncol = 1)
  pX <- prob_logistic(y = Si_theta[, X_val], mu = mu4)
  d_gamma7 <- (2 * Si_theta[, X_val] - 1) * sigmoid(mu4) * (1 - sigmoid(mu4))

  # Keep from dividing by 0
  pY[pY == 0] <- 1
  pYstar[pYstar == 0] <- 1
  pXstar[pXstar == 0] <- 1
  pX[pX == 0] <- 1

  # Validated
  Si_theta$Si_beta[val] <- d_beta[val] / pY[val]
  Si_theta$Si_alpha[val] <- d_alpha[val] / pY[val]
  Si_theta$Si_gamma0[val] <- d_gamma0[val] / pYstar[val]
  Si_theta$Si_gamma1[val] <- d_gamma1[val] / pYstar[val]
  Si_theta$Si_gamma2[val] <- d_gamma2[val] / pYstar[val]
  Si_theta$Si_gamma3[val] <- d_gamma3[val] / pYstar[val]
  Si_theta$Si_gamma4[val] <- d_gamma4[val] / pXstar[val]
  Si_theta$Si_gamma5[val] <- d_gamma5[val] / pXstar[val]
  Si_theta$Si_gamma6[val] <- d_gamma6[val] / pXstar[val]
  Si_theta$Si_gamma7[val] <- d_gamma7[val] / pX[val]

  ## Get the joint distributions of (Y, X, Y*, X*) -------------
  joint_v <- pY[val] * pYstar[val] * pXstar[val] * pX[val]
  joint_uv <- pY[-val] * pYstar[-val] * pXstar[-val] * pX[-val]

  # Save the individual scores of validated subjects
  Si_theta_v <- Si_theta[val, ]

  # Unvalidated
  Si_theta$PhI_strat <- paste0("(", Si_theta[, Y_unval], ",", Si_theta[, X_unval], ")")
  denom <- rowsum(x = joint_uv, group = Si_theta[-val, "PhI_strat"], reorder = F)

  ## (among unvalidated) Get the joint distributions of (Y*, X*) -
  phI_joint_uv <- rowsum(x = joint_uv, group = Si_theta[-val, "PhI_strat"], reorder = FALSE)

  Si_theta_uv <- Si_theta[-val, ]
  Si_theta_uv[, Y_val] <- Si_theta_uv[, X_val] <- NA
  Si_theta_uv <- unique(Si_theta_uv)

  # For the unvalidated subjects, the scores are only based on (Y*,X*)
  # So there should be just 4 rows returned with Validated = 0
  Si_theta_uv$Si_beta <- rowsum(x = d_beta[-val] * pYstar[-val] * pXstar[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_alpha <- rowsum(x = d_alpha[-val] * pYstar[-val] * pXstar[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_gamma0 <- rowsum(x = pY[-val] * d_gamma0[-val] * pXstar[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_gamma1 <- rowsum(x = pY[-val] * d_gamma1[-val] * pXstar[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_gamma2 <- rowsum(x = pY[-val] * d_gamma2[-val] * pXstar[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_gamma3 <- rowsum(x = pY[-val] * d_gamma3[-val] * pXstar[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_gamma4 <- rowsum(x = pY[-val] * pYstar[-val] * d_gamma4[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_gamma5 <- rowsum(x = pY[-val] * pYstar[-val] * d_gamma5[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_gamma6 <- rowsum(x = pY[-val] * pYstar[-val] * d_gamma6[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv$Si_gamma7 <- rowsum(x = pY[-val] * pYstar[-val] * pXstar[-val] * d_gamma7[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom

  Si_theta_uv$PhI_strat <- NULL

  Si_theta <- rbind(cbind(Si_theta_v, joint_exV = joint_v),
                    cbind(Si_theta_uv, joint_exV = phI_joint_uv))

  return(Si_theta)
}
