#' Sample unvalidated balanced case-control (BCC*) audit based on Phase I variables.
#' @name sample_bcc
#' @param dat Dataframe or matrix containing columns \code{Y_unval}, \code{Y_val}, \code{X_unval}, \code{X_val}.
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index). (Left \code{NULL} if \code{errors_in = "X only"}.)
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index). (Left \code{NULL} if \code{errors_in = "Y only"}.)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param errors_in Measurement error setting, options are \code{"Both"}, \code{"Y only"}, \code{"X only"}. Default is \code{"Both"}.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_bcc <- function(dat, Y_unval = NULL, Y_val, X_unval = NULL, X_val, phI, phII, errors_in = "Both") {

  if (errors_in == "Both") {
    # Cross-tabulate Phase I data (Y*,X*)
    phI_tab <- table(dat[, Y_unval], dat[, X_unval])

    # Create index vectors for subjects based on (Y*, X*)
    ind_N11 <- which(dat[, Y_unval] == 1 & dat[, X_unval] == 1)
    ind_N10 <- which(dat[, Y_unval] == 1 & dat[, X_unval] == 0)
    ind_N01 <- which(dat[, Y_unval] == 0 & dat[, X_unval] == 1)
    ind_N00 <- which(dat[, Y_unval] == 0 & dat[, X_unval] == 0)
  } else if (errors_in == "Y only") {
    # Cross-tabulate Phase I data (Y*,X)
    phI_tab <- table(dat[, Y_unval], dat[, X_val])

    # Create index vectors for subjects based on (Y*, X)
    ind_N11 <- which(dat[, Y_unval] == 1 & dat[, X_val] == 1)
    ind_N10 <- which(dat[, Y_unval] == 1 & dat[, X_val] == 0)
    ind_N01 <- which(dat[, Y_unval] == 0 & dat[, X_val] == 1)
    ind_N00 <- which(dat[, Y_unval] == 0 & dat[, X_val] == 0)
  } else if (errors_in == "X only") {
    # Cross-tabulate Phase I data (Y,X*)
    phI_tab <- table(dat[, Y_val], dat[, X_unval])

    # Create index vectors for subjects based on (Y, X*)
    ind_N11 <- which(dat[, Y_val] == 1 & dat[, X_unval] == 1)
    ind_N10 <- which(dat[, Y_val] == 1 & dat[, X_unval] == 0)
    ind_N01 <- which(dat[, Y_val] == 0 & dat[, X_unval] == 1)
    ind_N00 <- which(dat[, Y_val] == 0 & dat[, X_unval] == 0)
  }

  N00 <- phI_tab[1,1]; N01 <- phI_tab[1,2]; N10 <- phI_tab[2,1]; N11 <- phI_tab[2,2]
  phI_strat <- list(N00 = N00, N01 = N01, N10 = N10, N11 = N11)

  # If data are available, take n / 4 from each of the strata
  aud <- pmin(c(phI_strat$N11, phI_strat$N10, phI_strat$N01, phI_strat$N00), rep(phII / 4, 4))

  # If any strata were too small to take n / 4, allocate remaining audit evenly between remaining strata
  while (sum(aud) < phII) {
    if (aud[1] < phI_strat$N11) {
      aud[1] <- aud[1] + 1
      if (sum(aud) == phII) break()
    }
    if (aud[2] < phI_strat$N10) {
      aud[2] <- aud[2] + 1
      if (sum(aud) == phII) break()
    }
    if (aud[3] < phI_strat$N01) {
      aud[3] <- aud[3] + 1
      if (sum(aud) == phII) break()
    }
    if (aud[4] < phI_strat$N00) {
      aud[4] <- aud[4] + 1
      if (sum(aud) == phII) break()
    }
  }

  V_bcc <- as.numeric(seq(1, phI) %in% c(sample(x = ind_N11, size = aud[1], replace = F),
                                         sample(x = ind_N10, size = aud[2], replace = F),
                                         sample(x = ind_N01, size = aud[3], replace = F),
                                         sample(x = ind_N00, size = aud[4], replace = F)))
  return(V_bcc)
}
