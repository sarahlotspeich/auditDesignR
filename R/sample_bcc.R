#' Sample unvalidated balanced case-control (BCC*) audit based on Phase I variables.
#' @name sample_bcc
#' @param dat Dataframe or matrix containing columns \code{sample_on}.
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 binary variables can be accommodated.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_bcc <- function(dat, Y_unval = NULL, Y_val, X_unval = NULL, X_val, phI, phII, sample_on) {

  # Cross-tabulate Phase I data (sample_on)
  phI_tab <- table(dat[, sample_on])
  num_strat <- 2 ^ length(sample_on)

  if (num_strat == 2) {
    # Create index vectors for subjects
    ind_N1 <- which(dat[, sample_on[1]] == 1)
    ind_N0 <- which(dat[, sample_on[1]] == 0)

    # And stratum sizes for full cohort
    N1 <- length(ind_N1)
    N0 <- length(ind_N0)

    phI_strat <- list(N0 = N0, N1 = N1)
  } else if (num_strat == 4) {
    # Create index vectors for subjects
    ind_N11 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1)
    ind_N10 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0)
    ind_N01 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1)
    ind_N00 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0)

    # And stratum sizes for full cohort
    N11 <- length(ind_N11)
    N10 <- length(ind_N10)
    N01 <- length(ind_N01)
    N00 <- length(ind_N00)

    phI_strat <- list(N00 = N00, N01 = N01, N10 = N10, N11 = N11)
  } else if (num_strat == 8) {
    # Create index vectors for subjects
    ind_N11_0 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 0)
    ind_N10_0 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 0)
    ind_N01_0 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 0)
    ind_N00_0 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 0)
    ind_N11_1 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 1)
    ind_N10_1 <- which(dat[, sample_on[1]] == 1 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 1)
    ind_N01_1 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 1 & dat[, sample_on[3]] == 1)
    ind_N00_1 <- which(dat[, sample_on[1]] == 0 & dat[, sample_on[2]] == 0 & dat[, sample_on[3]] == 1)

    # And stratum sizes for full cohort
    N11_0 <- length(ind_N11_0)
    N10_0 <- length(ind_N10_0)
    N01_0 <- length(ind_N01_0)
    N00_0 <- length(ind_N00_0)
    N11_1 <- length(ind_N11_1)
    N10_1 <- length(ind_N10_1)
    N01_1 <- length(ind_N01_1)
    N00_1 <- length(ind_N00_1)

    phI_strat <- list(N00_0 = N00_0, N01_0 = N01_0, N10_0 = N10_0, N11_0 = N11_0,
                      N00_1 = N00_1, N01_1 = N01_1, N10_1 = N10_1, N11_1 = N11_1)
  }

  # If data are available, take phII / num_strat from each of the strata
  aud <- pmin(unlist(phI_strat), floor(phII / num_strat))

  # If any strata were too small to take phII / num_strat, allocate remaining audit evenly between remaining strata
  if (num_strat == 2) {
    while (sum(aud) < phII) {
      if (aud$N1 < phI_strat$N1) {
        aud$N1 <- aud$N1 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N0 < phI_strat$N0) {
        aud$N0 <- aud$N0 + 1
        if (sum(aud) == phII) break()
      }
    }

    V_bcc <- as.numeric(seq(1, phI) %in% c(sample(x = ind_N0, size = aud$N0, replace = F),
                                           sample(x = ind_N1, size = aud$N1, replace = F)))
  } else if (num_strat == 4) {
    while (sum(aud) < phII) {
      if (aud$N11 < phI_strat$N11) {
        aud$N11 <- aud$N11 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N10 < phI_strat$N10) {
        aud$N10 <- aud$N10 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N01 < phI_strat$N01) {
        aud$N01 <- aud$N01 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N00 < phI_strat$N00) {
        aud$N00 <- aud$N00 + 1
        if (sum(aud) == phII) break()
      }
    }

    V_bcc <- as.numeric(seq(1, phI) %in% c(sample(x = ind_N11, size = aud$N11, replace = F),
                                           sample(x = ind_N10, size = aud$N10, replace = F),
                                           sample(x = ind_N01, size = aud$N01, replace = F),
                                           sample(x = ind_N00, size = aud$N00, replace = F)))
  } else if (num_strat == 8) {
    while (sum(aud) < phII) {
      if (aud$N11_1 < phI_strat$N11_1) {
        aud$N11_1 <- aud$N11_1 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N11_0 < phI_strat$N11_0) {
        aud$N11_0 <- aud$N11_0 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N10_1 < phI_strat$N10_1) {
        aud$N10_1 <- aud$N10_1 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N10_0 < phI_strat$N10_0) {
        aud$N10_0 <- aud$N10_0 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N01_1 < phI_strat$N01_1) {
        aud$N01_1 <- aud$N01_1 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N01_0 < phI_strat$N01_0) {
        aud$N01_0 <- aud$N01_0 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N00_1 < phI_strat$N00_1) {
        aud$N00_1 <- aud$N00_1 + 1
        if (sum(aud) == phII) break()
      }
      if (aud$N00_0 < phI_strat$N00_0) {
        aud$N00_0 <- aud$N00_0 + 1
        if (sum(aud) == phII) break()
      }
    }

    V_bcc <- as.numeric(seq(1, phI) %in% c(sample(x = ind_N11_1, size = aud$N11_1, replace = F),
                                           sample(x = ind_N10_1, size = aud$N10_1, replace = F),
                                           sample(x = ind_N01_1, size = aud$N01_1, replace = F),
                                           sample(x = ind_N00_1, size = aud$N00_1, replace = F),
                                           sample(x = ind_N11_0, size = aud$N11_0, replace = F),
                                           sample(x = ind_N10_0, size = aud$N10_0, replace = F),
                                           sample(x = ind_N01_0, size = aud$N01_0, replace = F),
                                           sample(x = ind_N00_0, size = aud$N00_0, replace = F)))
  }
  return(V_bcc)
}
