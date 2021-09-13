score_w_int <- function(comp_dat, Y_val, Y_unval = NULL, X_val, X_unval = NULL, addl_covar = NULL, Validated, 
                        nondiff_Y_unval = FALSE, nondiff_X_unval = FALSE, beta, eta,
                        int_Y_val = FALSE, int_Y_unval = FALSE, int_X_val = FALSE, int_X_unval = FALSE) {
  val <- which(comp_dat[, Validated] == 1)
  n <- length(val)

  full_eta <- eta

  # Parameters P(Y_val|X_val, addl_covar) - coeff for X_val were already saved as beta
  if (int_Y_val) {
    alpha <- eta[1:(2 + length(addl_covar))]  
    eta <- eta[-c(1:(2 + length(addl_covar)))]
  } else {
    alpha <- eta[1:(1 + length(addl_covar))]  
    eta <- eta[-c(1:(1 + length(addl_covar)))]
  }

  # Parameters P(Y_unval|X_unval, Y_val, X_val, addl_covar)
  if (!is.null(Y_unval)) {
    if (!nondiff_Y_unval) {
      if (int_Y_unval) {
        gamma_Ystar <- eta[1:(2 + length(c(Y_val, X_val, X_unval, addl_covar)))]
        eta <- eta[-c(1:(2 + length(c(Y_val, X_val, X_unval, addl_covar))))]
      } else {
        gamma_Ystar <- eta[1:(1 + length(c(Y_val, X_val, X_unval, addl_covar)))]
        eta <- eta[-c(1:(1 + length(c(Y_val, X_val, X_unval, addl_covar))))]
      }
    } else {
      if (int_Y_unval) {
        gamma_Ystar <- eta[1:(2 + length(c(Y_val, addl_covar)))]
        eta <- eta[-c(1:(2 + length(c(Y_val, addl_covar))))]
      } else {
        gamma_Ystar <- eta[1:(1 + length(c(Y_val, addl_covar)))]
        eta <- eta[-c(1:(1 + length(c(Y_val, addl_covar))))]
      }
    }
  }

  # Parameters P(X_unval|Y_val, X_val, addl_covar)
  if (!is.null(X_unval)) {
    if (!nondiff_X_unval) {
      if (int_X_unval) {
        gamma_Xstar <- eta[1:(2 + length(c(Y_val, X_val, addl_covar)))]
        eta <- eta[-c(1:(2 + length(c(Y_val, X_val, addl_covar))))]
      } else {
        gamma_Xstar <- eta[1:(1 + length(c(Y_val, X_val, addl_covar)))]
        eta <- eta[-c(1:(1 + length(c(Y_val, X_val, addl_covar))))]
      }
    } else {
      if (int_X_unval) {
        gamma_Xstar <- eta[1:(2 + length(c(X_val, addl_covar)))]
        eta <- eta[-c(1:(2 + length(c(X_val, addl_covar))))]
      } else {
        gamma_Xstar <- eta[1:(1 + length(c(X_val, addl_covar)))]
        eta <- eta[-c(1:(1 + length(c(X_val, addl_covar))))]
      }
    }
  }

  # Parameters P(X_val|addl_covar)
  if (int_X_val) {
    gamma_X <- eta[1:(2 + length(addl_covar))]
  } else {
    gamma_X <- eta[1:(1 + length(addl_covar))]
  }

  Si_theta <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(c(beta, full_eta)))
  colnames(Si_theta) <- c("Si_beta", paste0("Si_gamma", seq(1, length(full_eta))))

  Si_theta <- data.frame(comp_dat, Si_theta, row.names = NULL)

  ## S(beta), S(alpha)
  if (int_Y_val) {
    XZ <- Si_theta[, addl_covar] * Si_theta[, X_val]
    mu1 <- data.matrix(cbind(1, Si_theta[, addl_covar], XZ, Si_theta[, X_val])) %*% matrix(c(alpha, beta), ncol = 1)
    pY <- prob_logistic(y = Si_theta[, Y_val], mu = mu1)
    d_beta <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(c(alpha, beta)))
    d_beta[, 1] <- (2 * Si_theta[, Y_val] - 1) * sigmoid(mu1) * (1 - sigmoid(mu1))
    t <- cbind(Si_theta[, addl_covar], XZ, Si_theta[, X_val])
    for (p in 2:ncol(d_beta)) {
      d_beta[, p] <- t[, (p - 1)] * d_beta[, 1]
    }
  } else {
    mu1 <- data.matrix(cbind(1, Si_theta[, c(addl_covar, X_val)])) %*% matrix(c(alpha, beta), ncol = 1)
    pY <- prob_logistic(y = Si_theta[, Y_val], mu = mu1)
    d_beta <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(c(alpha, beta)))
    d_beta[, 1] <- (2 * Si_theta[, Y_val] - 1) * sigmoid(mu1) * (1 - sigmoid(mu1))
    for (p in 2:ncol(d_beta)) {
      d_beta[, p] <- Si_theta[, c(addl_covar, X_val)[(p - 1)]] * d_beta[, 1]
    }
  }
  
  ## S(gamma0), ... , S(gamma3)
  if (!is.null(Y_unval)) {
    if (!nondiff_Y_unval) {
      if (int_Y_unval) {
        XZ <- Si_theta[, addl_covar] * Si_theta[, X_val]
        mu2 <- data.matrix(cbind(1, Si_theta[, c(X_unval, Y_val, X_val, addl_covar)], XZ)) %*% matrix(gamma_Ystar, ncol = 1)
        pYstar <- prob_logistic(y = Si_theta[, Y_unval], mu = mu2)
        d_gamma <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_Ystar))
        d_gamma[, 1] <- (2 * Si_theta[, Y_unval] - 1) * sigmoid(mu2) * (1 - sigmoid(mu2))
        t <- cbind(Si_theta[, c(X_unval, Y_val, X_val, addl_covar)], XZ)
        for (p in 2:ncol(d_gamma)) {
          d_gamma[, p] <- t[, (p - 1)] * d_gamma[, 1]
        }
      } else {
        mu2 <- data.matrix(cbind(1, Si_theta[, c(X_unval, Y_val, X_val, addl_covar)])) %*% matrix(gamma_Ystar, ncol = 1)
        pYstar <- prob_logistic(y = Si_theta[, Y_unval], mu = mu2)
        d_gamma <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_Ystar))
        d_gamma[, 1] <- (2 * Si_theta[, Y_unval] - 1) * sigmoid(mu2) * (1 - sigmoid(mu2))
        for (p in 2:ncol(d_gamma)) {
          d_gamma[, p] <- Si_theta[, c(X_unval, Y_val, X_val, addl_covar)[(p - 1)]] * d_gamma[, 1]
        }
      }
    } else {
      if (int_Y_unval) {
        XZ <- Si_theta[, addl_covar] * Si_theta[, X_val]
        mu2 <- data.matrix(cbind(1, Si_theta[, c(Y_val, addl_covar)], XZ)) %*% matrix(gamma_Ystar, ncol = 1)
        pYstar <- prob_logistic(y = Si_theta[, Y_unval], mu = mu2)
        d_gamma <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_Ystar))
        d_gamma[, 1] <- (2 * Si_theta[, Y_unval] - 1) * sigmoid(mu2) * (1 - sigmoid(mu2))
        t <- cbind(Si_theta[, c(Y_val, addl_covar)], XZ)
        for (p in 2:ncol(d_gamma)) {
          d_gamma[, p] <- t[, (p - 1)] * d_gamma[, 1]
        }
      } else {
        mu2 <- data.matrix(cbind(1, Si_theta[, c(Y_val, addl_covar)])) %*% matrix(gamma_Ystar, ncol = 1)
        pYstar <- prob_logistic(y = Si_theta[, Y_unval], mu = mu2)
        d_gamma <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_Ystar))
        d_gamma[, 1] <- (2 * Si_theta[, Y_unval] - 1) * sigmoid(mu2) * (1 - sigmoid(mu2))
        for (p in 2:ncol(d_gamma)) {
          d_gamma[, p] <- Si_theta[, c(Y_val, addl_covar)[(p - 1)]] * d_gamma[, 1]
        }
      }
    }
  } else { pYstar <- rep(1, nrow(Si_theta)) }

  ## S(gamma4), ..., S(gamma6)
  if (!is.null(X_unval)) {
    if (!nondiff_X_unval) {
      if (int_X_unval) {
        XZ <- Si_theta[, addl_covar] * Si_theta[, X_val]
        mu3 <- data.matrix(cbind(1, Si_theta[, c(Y_val, X_val, addl_covar)], XZ)) %*% matrix(gamma_Xstar, ncol = 1)
        pXstar <- prob_logistic(y = Si_theta[, X_unval], mu = mu3)
        d_delta <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_Xstar))
        d_delta[, 1] <- (2 * Si_theta[, X_unval] - 1) * sigmoid(mu3) * (1 - sigmoid(mu3))
        t <- cbind(Si_theta[, c(Y_val, X_val, addl_covar)], XZ)
        for (p in 2:ncol(d_delta)) {
          d_delta[, p] <- t[, (p - 1)] * d_delta[, 1]
        }
      } else {
        mu3 <- data.matrix(cbind(1, Si_theta[, c(Y_val, X_val, addl_covar)])) %*% matrix(gamma_Xstar, ncol = 1)
        pXstar <- prob_logistic(y = Si_theta[, X_unval], mu = mu3)
        d_delta <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_Xstar))
        d_delta[, 1] <- (2 * Si_theta[, X_unval] - 1) * sigmoid(mu3) * (1 - sigmoid(mu3))
        for (p in 2:ncol(d_delta)) {
          d_delta[, p] <- Si_theta[, c(Y_val, X_val, addl_covar)[(p - 1)]] * d_delta[, 1]
        }
      }
    } else {
      if (int_X_unval) {
        XZ <- Si_theta[, addl_covar] * Si_theta[, X_val]
        mu3 <- data.matrix(cbind(1, Si_theta[, c(X_val, addl_covar)], XZ)) %*% matrix(gamma_Xstar, ncol = 1)
        pXstar <- prob_logistic(y = Si_theta[, X_unval], mu = mu3)
        d_delta <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_Xstar))
        d_delta[, 1] <- (2 * Si_theta[, X_unval] - 1) * sigmoid(mu3) * (1 - sigmoid(mu3))
        t <- cbind(Si_theta[, c(X_val, addl_covar)], XZ)
        for (p in 2:ncol(d_delta)) {
          d_delta[, p] <- t[, (p - 1)] * d_delta[, 1]
        }
      } else {
        mu3 <- data.matrix(cbind(1, Si_theta[, c(X_val, addl_covar)])) %*% matrix(gamma_Xstar, ncol = 1)
        pXstar <- prob_logistic(y = Si_theta[, X_unval], mu = mu3)
        d_delta <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_Xstar))
        d_delta[, 1] <- (2 * Si_theta[, X_unval] - 1) * sigmoid(mu3) * (1 - sigmoid(mu3))
        for (p in 2:ncol(d_delta)) {
          d_delta[, p] <- Si_theta[, c(X_val, addl_covar)[(p - 1)]] * d_delta[, 1]
        }
      }
    }
  }  else { pXstar <- rep(1, nrow(Si_theta)) }

  ## S(gamma7)
  if (!is.null(addl_covar)) {
    if (int_X_val) {
      mu4 <- data.matrix(cbind(1, Si_theta[, c(addl_covar)], XZ)) %*% matrix(gamma_X, ncol = 1)
    } else {
      mu4 <- data.matrix(cbind(1, Si_theta[, c(addl_covar)])) %*% matrix(gamma_X, ncol = 1) 
    }
  } else {
    mu4 <- matrix(rep(gamma_X, nrow(Si_theta)), ncol = 1)
  }
  pX <- prob_logistic(y = Si_theta[, X_val], mu = mu4)
  d_tau <- matrix(data = NA, nrow = nrow(comp_dat), ncol = length(gamma_X))
  d_tau[, 1] <- (2 * Si_theta[, X_val] - 1) * sigmoid(mu4) * (1 - sigmoid(mu4))
  if (int_X_val) {
    t <- cbind(Si_theta[, c(addl_covar)], XZ)
  } else {
    t <- cbind(Si_theta[, c(addl_covar)])
  }
  if (ncol(d_tau) > 1) {
    for (p in 2:ncol(d_tau)) {
      d_tau[, p] <- t[, (p - 1)] * d_tau[, 1]
    }
  }

  # Keep from dividing by 0
  pY[pY == 0] <- 1
  pYstar[pYstar == 0] <- 1
  pXstar[pXstar == 0] <- 1
  pX[pX == 0] <- 1

  score_cols <- grep(pattern = "Si_", colnames(Si_theta)) # indices for score columns in Si_theta

  # Validated
  Si_theta[val, score_cols[1]] <- (d_beta[val, ] / pY[val])[, ncol(d_beta)] # d_beta[val] / pY[val]
  Si_theta[val, score_cols[2:length(c(alpha, beta))]] <- (d_beta[val, ] / pY[val])[, -ncol(d_beta)]
  score_cols <- score_cols[-c(1:length(c(alpha, beta)))]

  if (!is.null(Y_unval)) {
    Si_theta[val, score_cols[1:length(gamma_Ystar)]] <- (d_gamma[val, ] / pYstar[val])
    score_cols <- score_cols[-c(1:length(gamma_Ystar))]
  }

  if (!is.null(X_unval)) {
    Si_theta[val, score_cols[1:length(gamma_Xstar)]] <- (d_delta[val, ] / pXstar[val])
    score_cols <- score_cols[-c(1:length(gamma_Xstar))]
  }

  Si_theta[val, score_cols[1:length(gamma_X)]] <- (d_tau[val, ] / pX[val])

  ## Get the joint distributions of (Y, X, Y*, X*) -------------
  joint_v <- pY[val] * pYstar[val] * pXstar[val] * pX[val]
  joint_uv <- pY[-val] * pYstar[-val] * pXstar[-val] * pX[-val]

  # Save the individual scores of validated subjects
  Si_theta_v <- Si_theta[val, ]

  # Unvalidated
  if (!is.null(Y_unval) & !is.null(X_unval)) {
    PhI_vars <- c(Y_unval, X_unval, addl_covar)
  } else if (!is.null(Y_unval)) {
    PhI_vars <- c(Y_unval, X_val, addl_covar)
  } else if (!is.null(X_unval)) {
    PhI_vars <- c(Y_val, X_unval, addl_covar)
  }
  Si_theta$PhI_strat <- apply(X = Si_theta[, PhI_vars], MARGIN = 1, FUN = paste, collapse = ",")
  denom <- rowsum(x = joint_uv, group = Si_theta[-val, "PhI_strat"], reorder = F)

  ## (among unvalidated) Get the joint distributions of (Y*, X*, Z) -
  phI_joint_uv <- rowsum(x = joint_uv, group = Si_theta[-val, "PhI_strat"], reorder = FALSE)

  Si_theta_uv <- Si_theta[-val, ]
  if (!is.null(Y_unval)) { Si_theta_uv[, Y_val] <- NA }
  if (!is.null(X_unval)) { Si_theta_uv[, X_val] <- NA }
  Si_theta_uv <- unique(Si_theta_uv)

  score_cols <- grep(pattern = "Si_", colnames(Si_theta_uv)) # indices for score columns in Si_theta

  # For the unvalidated subjects, the scores are only based on (Y*,X*,Z)
  # So there should be just 4/8 rows returned with Validated = 0

  Si_theta_uv[, score_cols[1]] <- rowsum(x = d_beta[-val, ncol(d_beta)] * pYstar[-val] * pXstar[-val] * pX[-val], group = Si_theta[-val, "PhI_strat"], reorder = F) / denom
  Si_theta_uv[, score_cols[2:length(c(alpha, beta))]] <- rowsum(x = d_beta[-val, -ncol(d_beta)] * c(pYstar[-val] * pXstar[-val] * pX[-val]), group = Si_theta[-val, "PhI_strat"], reorder = F) / c(denom)
  score_cols <- score_cols[-c(1:length(c(alpha, beta)))]

  if (!is.null(Y_unval)) {
    Si_theta_uv[, score_cols[1:length(gamma_Ystar)]] <- rowsum(x = d_gamma[-val, ] * c(pY[-val] * pXstar[-val] * pX[-val]), group = Si_theta[-val, "PhI_strat"], reorder = F) / c(denom)
    score_cols <- score_cols[-c(1:length(gamma_Ystar))]
  }

  if (!is.null(X_unval)) {
    Si_theta_uv[, score_cols[1:length(gamma_Xstar)]] <- rowsum(x = d_delta[-val, ] * c(pY[-val] * pYstar[-val] * pX[-val]), group = Si_theta[-val, "PhI_strat"], reorder = F) / c(denom)
    score_cols <- score_cols[-c(1:length(gamma_Xstar))]
  }

  Si_theta_uv[, score_cols[1:length(gamma_X)]] <- rowsum(x = c(pY[-val] * pYstar[-val] * pXstar[-val]) * d_tau[-val, ], group = Si_theta[-val, "PhI_strat"], reorder = F) / c(denom)

  Si_theta_uv$PhI_strat <- NULL

  Si_theta <- rbind(cbind(Si_theta_v, joint_exV = joint_v),
                    cbind(Si_theta_uv, joint_exV = phI_joint_uv))

  return(Si_theta)
}
