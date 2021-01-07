#' Sample optimal MLE (optMLE) audit based on Phase I variables.
#' @name sample_optMLE
#' @param dat Dataframe or matrix containing columns \code{Y_unval}, \code{Y_val}, \code{X_unval}, \code{X_val}.
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index). (Left \code{NULL} if \code{errors_in = "X only"}.)
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index). (Left \code{NULL} if \code{errors_in = "Y only"}.)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param phI Phase I sample size.
#' @param des Dataframe containing the 4 stratum sample sizes.
#' @param errors_in Measurement error setting, options are \code{"Both"}, \code{"Y only"}, \code{"X only"}. Default is \code{"Both"}.
#' @param wave1_Validated (For use with multi-wave designs.) Columns with the validation indicator (can be name or numeric index) from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_optMLE <- function(dat, Y_unval = NULL, Y_val, X_unval = NULL, X_val, phI, des, errors_in = "Both", wave1_Validated = NULL) {
  if (errors_in == "Both") {
    if (!is.null(wave1_Validated)) {
      # Create index vectors for subjects based on (Y*, X*)
      ind_N11 <- which(dat[, Y_unval] == 1 & dat[, X_unval] == 1 & dat[, wave1_Validated] == 0)
      ind_N10 <- which(dat[, Y_unval] == 1 & dat[, X_unval] == 0 & dat[, wave1_Validated] == 0)
      ind_N01 <- which(dat[, Y_unval] == 0 & dat[, X_unval] == 1 & dat[, wave1_Validated] == 0)
      ind_N00 <- which(dat[, Y_unval] == 0 & dat[, X_unval] == 0 & dat[, wave1_Validated] == 0)
    } else {
      # Create index vectors for subjects based on (Y*, X*)
      ind_N11 <- which(dat[, Y_unval] == 1 & dat[, X_unval] == 1)
      ind_N10 <- which(dat[, Y_unval] == 1 & dat[, X_unval] == 0)
      ind_N01 <- which(dat[, Y_unval] == 0 & dat[, X_unval] == 1)
      ind_N00 <- which(dat[, Y_unval] == 0 & dat[, X_unval] == 0)
    }
  } else if (errors_in == "Y only") {
    if (!is.null(wave1_Validated)) {
      # Create index vectors for subjects based on (Y*, X)
      ind_N11 <- which(dat[, Y_unval] == 1 & dat[, X_val] == 1 & dat[, wave1_Validated] == 0)
      ind_N10 <- which(dat[, Y_unval] == 1 & dat[, X_val] == 0 & dat[, wave1_Validated] == 0)
      ind_N01 <- which(dat[, Y_unval] == 0 & dat[, X_val] == 1 & dat[, wave1_Validated] == 0)
      ind_N00 <- which(dat[, Y_unval] == 0 & dat[, X_val] == 0 & dat[, wave1_Validated] == 0)
    } else {
      # Create index vectors for subjects based on (Y*, X)
      ind_N11 <- which(dat[, Y_unval] == 1 & dat[, X_val] == 1)
      ind_N10 <- which(dat[, Y_unval] == 1 & dat[, X_val] == 0)
      ind_N01 <- which(dat[, Y_unval] == 0 & dat[, X_val] == 1)
      ind_N00 <- which(dat[, Y_unval] == 0 & dat[, X_val] == 0)
    }
  } else if (errors_in == "X only") {
    if (!is.null(wave1_Validated)) {
      # Create index vectors for subjects based on (Y, X*)
      ind_N11 <- which(dat[, Y_val] == 1 & dat[, X_unval] == 1 & dat[, wave1_Validated] == 0)
      ind_N10 <- which(dat[, Y_val] == 1 & dat[, X_unval] == 0 & dat[, wave1_Validated] == 0)
      ind_N01 <- which(dat[, Y_val] == 0 & dat[, X_unval] == 1 & dat[, wave1_Validated] == 0)
      ind_N00 <- which(dat[, Y_val] == 0 & dat[, X_unval] == 0 & dat[, wave1_Validated] == 0)
    } else {
      # Create index vectors for subjects based on (Y, X*)
      ind_N11 <- which(dat[, Y_val] == 1 & dat[, X_unval] == 1)
      ind_N10 <- which(dat[, Y_val] == 1 & dat[, X_unval] == 0)
      ind_N01 <- which(dat[, Y_val] == 0 & dat[, X_unval] == 1)
      ind_N00 <- which(dat[, Y_val] == 0 & dat[, X_unval] == 0)
    }
  }
  V_optMLE <- as.numeric(seq(1, phI) %in% c(sample(x = ind_N11, size = des$n11, replace = F),
                                            sample(x = ind_N10, size = des$n10, replace = F),
                                            sample(x = ind_N01, size = des$n01, replace = F),
                                            sample(x = ind_N00, size = des$n00, replace = F)))
  return(V_optMLE)
}
