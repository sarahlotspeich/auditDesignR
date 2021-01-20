#' Information matrix for two-phase measurement error problems.
#' @name info
#' @param pi_probs Dataframe of sampling probabilities, including \code{sample_on}.
#' @param indiv_score Matrix of score vectors for all parameters, \code{beta} and \code{eta}.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index).
#' @return Information matrix.
#' @export
info <- function(pi_probs, indiv_score, sample_on) {

  # Based on pi_probs calculate P(V) for all Y*/X*
  pi_probs_full <- rbind(cbind(pi_probs[, sample_on], V = 1, pV = pi_probs[, "pV1"]),
                         cbind(pi_probs[, sample_on], V = 0, pV = (1 - pi_probs[, "pV1"])))

  # Merge P(V) into the indiv_score
  indiv_score_full <- merge(x = indiv_score, y = pi_probs_full)

  # Calculate the full joint density
  ## P(Y*,X*,Y,X,V) = P(V=1|Y*,X*)P(Y*,X*,Y,X) for validated subjects
  ## P(Y*,X*,V) = P(V=0|Y*,X*)P(Y*,X*) for unvalidated subjects
  indiv_score_full$joint_inclV <- indiv_score_full$pV * indiv_score_full$joint_exV

  # Subset to the score functino columns of indiv_score
  S_theta <- indiv_score_full[, grep(pattern = "Si_", colnames(indiv_score_full))]

  # Create empty matrix I_theta for the information functions
  I_theta <- matrix(0, nrow = ncol(S_theta), ncol = ncol(S_theta))

  # For each column in S_theta (each parameter)
  for (r in 1:nrow(I_theta)) {
    for (c in r:ncol(I_theta)) {
      I_theta[r, c] <- I_theta[c, r] <- sum(S_theta[, r] * S_theta[, c] * indiv_score_full$joint_inclV)
    }
  }

  return(I_theta)
}
