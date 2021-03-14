#' Derived asymptotic variance of the MLE for beta, the logOR of X on Y
#' @name var_formula
#' @param pi_vec Vector of sampling probabilities (returned from \code{calc_pi_vec}).
#' @param phI Phase I sample size.
#' @param indiv_score Matrix of score vectors for all parameters, \code{beta} and \code{eta} (returned from \code{score}).
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 variables can be accommodated.
#' @return Scalar function value for Var(beta).
#' @export
var_formula <- function(pi_vec, phI, indiv_score, sample_on) {
  pi_probs <- unique(indiv_score[, c(sample_on, "strat")])
  pi_probs <-  cbind(pi_probs, pV1 = unlist(c(pi_vec[gsub("N", "pi", pi_probs$strat)])))

  i <- info(indiv_score = indiv_score, pi_probs = pi_probs, sample_on = sample_on)
  cov <- tryCatch(expr = solve(i[1, 1] - i[1, -1] %*% solve(i[-1, -1]) %*% i[-1, 1]), error = function(err) {matrix(NA, nrow = nrow(i), ncol = ncol(i))})
  v_beta <- diag(cov)[1] / phI
  return(v_beta)
}
