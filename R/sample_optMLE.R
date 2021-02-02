#' Sample optimal MLE (optMLE) audit based on Phase I variables.
#' @name sample_optMLE
#' @param dat Dataframe or matrix containing columns \code{sample_on}.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 binary variables can be accommodated.
#' @param des Dataframe containing the stratum sample sizes.
#' @param wave1_Validated (For use with multi-wave designs.) Columns with the validation indicator (can be name or numeric index) from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_optMLE <- function(dat, des, errors_in = "Both", wave1_Validated = NULL) {
  phI <- nrow(dat)

  # Cross-tabulate Phase I data (sample_on)
  phI_tab <- table(dat[, sample_on])
  num_strat <- 2 ^ length(sample_on)

  if (num_strat == 2) {
    if (!is.null(wave1_Validated)) {
      # Create index vectors for subjects
      ## Exclude subjects who were already validated in wave 1
      ind_N1 <- which(dat[, sample_on[1]] == 1 & dat[, wave1_Validated] == 0)
      ind_N0 <- which(dat[, sample_on[1]] == 0 & dat[, wave1_Validated] == 0)
    } else {
      # Create index vectors for subjects
      ind_N1 <- which(dat[, sample_on[1]] == 1)
      ind_N0 <- which(dat[, sample_on[1]] == 0)
    }

    V_optMLE <- as.numeric(seq(1, phI) %in% c(sample(x = ind_N1, size = des$n1, replace = F),
                                              sample(x = ind_N0, size = des$n0, replace = F)))
  } else if (num_strat == 4) {
    if (!is.null(wave1_Validated)) {
      # Create index vectors for subjects
      ## Exclude subjects who were already validated in wave 1
      ind_N11 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1 & dat[, wave1_Validated] == 0)
      ind_N10 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0 & dat[, wave1_Validated] == 0)
      ind_N01 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1 & dat[, wave1_Validated] == 0)
      ind_N00 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0 & dat[, wave1_Validated] == 0)
    } else {
      # Create index vectors for subjects
      ind_N11 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1)
      ind_N10 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0)
      ind_N01 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1)
      ind_N00 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0)
    }

    V_optMLE <- as.numeric(seq(1, phI) %in% c(sample(x = ind_N11, size = des$n11, replace = F),
                                              sample(x = ind_N10, size = des$n10, replace = F),
                                              sample(x = ind_N01, size = des$n01, replace = F),
                                              sample(x = ind_N00, size = des$n00, replace = F)))
  } else if (num_strat == 8) {
    if (!is.null(wave1_Validated)) {
      # Create index vectors for subjects
      ## Exclude subjects who were already validated in wave 1
      ind_N11_0 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 0 & dat[, wave1_Validated] == 0)
      ind_N10_0 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 0 & dat[, wave1_Validated] == 0)
      ind_N01_0 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 0 & dat[, wave1_Validated] == 0)
      ind_N00_0 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 0 & dat[, wave1_Validated] == 0)
      ind_N11_1 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 1 & dat[, wave1_Validated] == 0)
      ind_N10_1 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 1 & dat[, wave1_Validated] == 0)
      ind_N01_1 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 1 & dat[, wave1_Validated] == 0)
      ind_N00_1 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 1 & dat[, wave1_Validated] == 0)
    } else {
      # Create index vectors for subjects
      ind_N11_0 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 0)
      ind_N10_0 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 0)
      ind_N01_0 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 0)
      ind_N00_0 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 0)
      ind_N11_1 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 1)
      ind_N10_1 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 1)
      ind_N01_1 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 1)
      ind_N00_1 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 1)
    }
    V_optMLE <- as.numeric(seq(1, phI) %in% c(sample(x = ind_N11_0, size = des$n11_0, replace = F),
                                              sample(x = ind_N10_0, size = des$n10_0, replace = F),
                                              sample(x = ind_N01_0, size = des$n01_0, replace = F),
                                              sample(x = ind_N00_0, size = des$n00_0, replace = F),
                                              sample(x = ind_N11_1, size = des$n11_1, replace = F),
                                              sample(x = ind_N10_1, size = des$n10_1, replace = F),
                                              sample(x = ind_N01_1, size = des$n01_1, replace = F),
                                              sample(x = ind_N00_1, size = des$n00_1, replace = F)))

  }
  return(V_optMLE)
}
