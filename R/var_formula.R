#' Derived asymptotic variance of the MLE for beta, the logOR of X on Y
#' @name var_formula
#' @param pi_vec Vector of sampling probabilities (returned from \code{calc_pi_vec}).
#' @param phI Phase I sample size.
#' @param indiv_score Matrix of score vectors for all parameters, \code{beta} and \code{eta} (returned from \code{score}).
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, \code{var_formula} can handle sampling on up to 3 variables.
#' @return Scalar function value for Var(beta).
#' @export
var_formula <- function(pi_vec = NULL, phI, indiv_score, sample_on) {
  if (length(sample_on) == 1) {
    pi_probs <- data.frame(Var1 = c(0, 1),
                           pV1 = c(pi_vec[1], pi_vec[2]))
    colnames(pi_probs)[1] <- sample_on
  } else if (length(sample_on) == 2) {
    pi_probs <- data.frame(Var1 = c(0, 0, 1, 1),
                           Var2 = c(0, 1, 0, 1),
                           pV1 = c(pi_vec[1], pi_vec[2], pi_vec[3], pi_vec[4]))
    colnames(pi_probs)[1:2] <- sample_on
  } else if (length(sample_on) == 3) {
    pi_probs <- data.frame(Var1 = c(0, 0, 1, 1, 0, 0, 1, 1),
                           Var2 = c(0, 1, 0, 1, 0, 1, 0, 1),
                           Var3 = c(0, 0, 0, 0, 1, 1, 1, 1),
                           pV1 = c(pi_vec[1], pi_vec[2], pi_vec[3], pi_vec[4],
                                   pi_vec[5], pi_vec[6], pi_vec[7], pi_vec[8]))
    colnames(pi_probs)[1:3] <- sample_on
  } else {
    return(warning("You are attempting to sample on too many variables."))
  }
  #
  # if (errors_in == "Both") {
  #   phI_vars <- c(Y_unval, X_unval)
  # } else if (errors_in == "X only") {
  #   phI_vars <- c(Y_val, X_unval)
  # } else if (errors_in == "Y only") {
  #   phI_vars <- c(Y_unval, X_val)
  # }
  # colnames(pi_probs)[1:2] <- phI_vars

  i <- info(indiv_score = indiv_score, pi_probs = pi_probs, sample_on = sample_on)
  cov <- tryCatch(expr = solve(i[1, 1] - i[1, -1] %*% solve(i[-1, -1]) %*% i[-1, 1]), error = function(err) {matrix(NA, nrow = nrow(i), ncol = ncol(i))})
  #cov <- tryCatch(expr = solve(i), error = function(err) {matrix(NA, nrow = nrow(i), ncol = ncol(i))})
  v_beta <- diag(cov)[1] / phI # round(diag(cov)[1] / phI, 8)
  #se_beta <- sqrt(v_beta)
  return(v_beta)
}
