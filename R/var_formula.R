#' Derived asymptotic variance of the MLE for beta, the logOR of X on Y
#' @name var_formula
#' @param pi_vec Vector of sampling probabilities (returned from \code{calc_pi_vec}).
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index). (Left \code{NULL} if \code{errors_in = "X only"}.)
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index). (Left \code{NULL} if \code{errors_in = "Y only"}.)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param addl_covar Column(s) with additional error-free covariates (can be name or numeric index).
#' @param indiv_score Matrix of score vectors for all parameters, \code{beta} and \code{eta}.
#' @param phI Phase I sample size.
#' @param errors_in Measurement error setting, options are \code{"Both"}, \code{"Y only"}, \code{"X only"}. Default is \code{"Both"}.
#' @return Scalar function value for Var(beta).
#' @export
var_formula <- function(pi_vec = NULL, Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, addl_covar = NULL, phI, indiv_score, errors_in = "Both") {

  pi_probs <- data.frame(Var1 = c(0, 0, 1, 1),
                         Var2 = c(0, 1, 0, 1),
                         pV1 = c(pi_vec[1], pi_vec[2], pi_vec[3], pi_vec[4]))

  if (errors_in == "Both") {
    phI_vars <- c(Y_unval, X_unval)
  } else if (errors_in == "X only") {
    phI_vars <- c(Y_val, X_unval)
  } else if (errors_in == "Y only") {
    phI_vars <- c(Y_unval, X_val)
  }
  colnames(pi_probs)[1:2] <- phI_vars

  i <- info(indiv_score = indiv_score, pi_probs = pi_probs, phI_vars = phI_vars)
  cov <- tryCatch(expr = solve(i[1, 1] - i[1, -1] %*% solve(i[-1, -1]) %*% i[-1, 1]), error = function(err) {matrix(NA, nrow = nrow(i), ncol = ncol(i))})
  #cov <- tryCatch(expr = solve(i), error = function(err) {matrix(NA, nrow = nrow(i), ncol = ncol(i))})
  v_beta <- diag(cov)[1] / phI # round(diag(cov)[1] / phI, 8)
  #se_beta <- sqrt(v_beta)
  return(v_beta)
}
