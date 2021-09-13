od_loglik_w_int <- function(params, dat, Y_val, Y_unval, X_val, X_unval, addl_covar, interact, Validated, 
                              int_Y_val = FALSE, int_Y_unval = FALSE, int_X_val = FALSE, int_X_unval = FALSE, for_nlm = TRUE) {
  N <- nrow(dat)
  dat[, "id"] <- seq(1, N)
  val <- which(as.numeric(dat[, Validated]) == 1)
  
  beta <- params[1:length(X_val)]
  eta <- params[-c(1:(length(X_val)))]
  
  if (int_Y_val) {
    # Parameters P(Y_val|X_val, addl_covar, int) - coeff for X_val were already saved as beta
    alpha <- eta[1:(1 + length(c(addl_covar, interact)))]
    eta <- eta[-c(1:(1 + length(c(addl_covar, interact))))]
    mu1 <- data.matrix(cbind(1, dat[val, c(addl_covar, interact, X_val)])) %*% matrix(c(alpha, beta), ncol = 1)
  } else {
    # Parameters P(Y_val|X_val, addl_covar) - coeff for X_val were already saved as beta
    alpha <- eta[1:(1 + length(addl_covar))]
    eta <- eta[-c(1:(1 + length(addl_covar)))]
    mu1 <- data.matrix(cbind(1, dat[val, c(addl_covar, X_val)])) %*% matrix(c(alpha, beta), ncol = 1)
  }
  
  if (int_Y_unval) {
    # Parameters P(Y_unval|X_unval, Y_val, X_val, addl_covar, int)
    gamma_Ystar <- eta[1:(1 + length(c(Y_val, X_val, X_unval, addl_covar, interact)))]
    eta <- eta[-c(1:(1 + length(c(Y_val, X_val, X_unval, addl_covar, interact))))]
    mu2 <- data.matrix(cbind(1, dat[val, c(X_unval, Y_val, X_val, addl_covar, interact)])) %*% matrix(gamma_Ystar, ncol = 1)
  } else {
    # Parameters P(Y_unval|X_unval, Y_val, X_val, addl_covar)
    gamma_Ystar <- eta[1:(1 + length(c(Y_val, X_val, X_unval, addl_covar)))]
    eta <- eta[-c(1:(1 + length(c(Y_val, X_val, X_unval, addl_covar))))]
    mu2 <- data.matrix(cbind(1, dat[val, c(X_unval, Y_val, X_val, addl_covar)])) %*% matrix(gamma_Ystar, ncol = 1)
  }
  
  if (int_X_unval) {
    # Parameters P(X_unval|Y_val, X_val, addl_covar, interact)
    gamma_Xstar <- eta[1:(1 + length(c(Y_val, X_val, addl_covar, interact)))]
    eta <- eta[-c(1:(1 + length(c(Y_val, X_val, addl_covar, interact))))]
    mu3 <- data.matrix(cbind(1, dat[val, c(Y_val, X_val, addl_covar, interact)])) %*% matrix(gamma_Xstar, ncol = 1)
  } else {
    # Parameters P(X_unval|Y_val, X_val, addl_covar)
    gamma_Xstar <- eta[1:(1 + length(c(Y_val, X_val, addl_covar)))]
    eta <- eta[-c(1:(1 + length(c(Y_val, X_val, addl_covar))))]
    mu3 <- data.matrix(cbind(1, dat[val, c(Y_val, X_val, addl_covar)])) %*% matrix(gamma_Xstar, ncol = 1)
  }
  
  if (int_X_val) {
    # Parameters P(X_val|addl_covar, interact)
    gamma_X <- eta[1:(1 + length(c(addl_covar, interact)))]
    mu4 <- data.matrix(cbind(1, dat[val, c(addl_covar, interact)])) %*% matrix(gamma_X, ncol = 1)
  } else {
    # Parameters P(X_val|addl_covar)
    gamma_X <- eta[1:(1 + length(addl_covar))]
    mu4 <- data.matrix(cbind(1, dat[val, c(addl_covar)])) %*% matrix(gamma_X, ncol = 1)
  }
  
  # Validated subjects ----------------------------------------------------------------------
  n <- length(val)
  
  pY <- prob_logistic(y = dat[val, Y_val], mu = mu1)
  pYstar <- prob_logistic(y = dat[val, Y_unval], mu = mu2)
  pXstar <- prob_logistic(y = dat[val, X_unval], mu = mu3)
  pX <- prob_logistic(y = dat[val, X_val], mu = mu4)
  
  l <- sum(log(pYstar) + log(pXstar) + log(pY) + log(pX))
  
  # Unvalidated subjects ---------------------------------------------------------------------
  unval <- which(as.numeric(dat[, Validated]) == 0)
  
  # Create complete dataset ------------------------------------------------------------------
  if (!is.null(X_unval)) {
    X_val_unique <- data.frame(unique(dat[val, X_val]))
    m <- nrow(X_val_unique)
    cd_unval <- dat[rep(unval, each = m), ]
    cd_unval[, X_val] <- X_val_unique[rep(seq(1, m), times = (N - n)), ]
    if (any(int_Y_val, int_Y_unval, int_X_unval, int_X_val)) {
      # Recalculate the interaction X * Z
      cd_unval[, interact] <- cd_unval[, X_val] * cd_unval[, addl_covar] 
    }
  } else {
    cd_unval <- dat[unval, ]
    m <- 1
  }
  
  if (!is.null(Y_unval)) {
    cd_unval <- rbind(cd_unval, cd_unval)
    cd_unval[, Y_val] <- rep(c(0, 1), each = ((N - n) * m))
  }
  
  if (int_Y_val) {
    mu1 <- data.matrix(cbind(1, cd_unval[, c(addl_covar, interact, X_val)])) %*% matrix(c(alpha, beta), ncol = 1) 
  } else {
    mu1 <- data.matrix(cbind(1, cd_unval[, c(addl_covar, X_val)])) %*% matrix(c(alpha, beta), ncol = 1)
  }
  pY <- prob_logistic(y = cd_unval[, Y_val], mu = mu1)
  
  if (int_Y_unval) {
    mu2 <- data.matrix(cbind(1, cd_unval[, c(X_unval, Y_val, X_val, addl_covar, interact)])) %*% matrix(gamma_Ystar, ncol = 1)
  } else {
    mu2 <- data.matrix(cbind(1, cd_unval[, c(X_unval, Y_val, X_val, addl_covar)])) %*% matrix(gamma_Ystar, ncol = 1)
  }
  pYstar <- prob_logistic(y = cd_unval[, Y_unval], mu = mu2)
  
  if (int_X_unval) {
    mu3 <- data.matrix(cbind(1, cd_unval[, c(Y_val, X_val, addl_covar, interact)])) %*% matrix(gamma_Xstar, ncol = 1)
  } else {
    mu3 <- data.matrix(cbind(1, cd_unval[, c(Y_val, X_val, addl_covar)])) %*% matrix(gamma_Xstar, ncol = 1)
  }
  pXstar <- prob_logistic(y = cd_unval[, X_unval], mu = mu3)
  
  if (int_X_val) {
    mu4 <- data.matrix(cbind(1, cd_unval[, c(addl_covar, interact)])) %*% matrix(gamma_X, ncol = 1)
  } else {
    mu4 <- data.matrix(cbind(1, cd_unval[, c(addl_covar)])) %*% matrix(gamma_X, ncol = 1)
  }
  pX <- prob_logistic(y = cd_unval[, X_val], mu = mu4)
  
  l <- l + sum(log(rowsum(x = pY * pYstar * pXstar * pX, group = cd_unval[, "id"])))
  
  if (for_nlm) { return(- l) } else { return(l) }
}
