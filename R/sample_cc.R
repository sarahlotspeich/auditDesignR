#' Sample unvalidated case-control (CC*) audit based on Phase I variables.
#' @name sample_cc
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
sample_cc <- function(dat, Y_unval = NULL, Y_val, X_unval = NULL, X_val, phI, phII, errors_in = "Both") {
  if (errors_in == "Both" || errors_in == "Y only") {
    # Cross-tabulate Phase I data (Y*)
    phI_tab <- table(dat[, Y_unval])

    cases <- phI_tab[1]
    controls <- phI_tab[2]

    # Create index vectors for subjects based on (Y*)
    ind_N0. <- which(dat[, Y_unval] == 0)
    ind_N1. <- which(dat[, Y_unval] == 1)
  } else if (errors_in == "X only") {
    # Cross-tabulate Phase I data (Y*)
    phI_tab <- table(dat[, Y_val])

    cases <- phI_tab[1]
    controls <- phI_tab[2]

    # Create index vectors for subjects based on (Y)
    ind_N0. <- which(dat[, Y_val] == 0)
    ind_N1. <- which(dat[, Y_val] == 1)
  }

  target_cases <- min(phII / 2, cases)
  target_controls <- min(phII / 2, controls)

  if (target_cases < phII / 2) {
    target_controls <- target_controls + (phII / 2 - target_cases)
  } else if (target_controls < phII / 2) {
    target_cases <- target_cases + (phII / 2 - target_controls)
  }

  # Randomly sample target_controls/target_cases from ind_N0. and ind_N1., respectively
  V_cc <- c(sample(x = ind_N1., size = target_cases, replace = F),
            sample(x = ind_N0., size = target_controls, replace = F))
  V_cc <- as.numeric(seq(1, phI) %in% V_cc[order(V_cc)])
  return(V_cc)
}
