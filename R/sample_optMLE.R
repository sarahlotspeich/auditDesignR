#' Sample optimal MLE (optMLE) audit based on Phase I variables.
#' @name sample_optMLE
#' @param dat Dataframe or matrix containing columns \code{sample_on}.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 binary variables can be accommodated.
#' @param des Dataframe containing the stratum sample sizes.
#' @param wave1_Validated (For use with multi-wave designs.) Columns with the validation indicator (can be name or numeric index) from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_optMLE <- function(dat, sample_on, des, wave1_Validated = NULL) {
  phI <- nrow(dat)
  # Cross-tabulate Phase I data (Y*)
  phI_tab <- data.frame(table(dat[, c(sample_on, wave1_Validated)]))
  colnames(phI_tab)[1:length(sample_on)] <- sample_on
  colnames(phI_tab)[ncol(phI_tab)] <- "N"
  #phI_tab$target <- pmin(floor(phII / nrow(phI_tab)), phI_tab$N)

  if (!is.null(wave1_Validated)) {
    phI_tab <- phI_tab[phI_tab[, wave1_Validated] == 0, ]
  }

  des_long <- data.frame(t(des[grep("n", colnames(des))]))
  colnames(des_long) <- "target"
  des_long$strat <- gsub("n", "N", rownames(des_long))

  # Create stratum IDs to sample on
  phI_tab$strat <- "N"
  dat$strat <- "N"
  for (s in sample_on) {
    phI_tab$strat <- paste0(phI_tab$strat, phI_tab[, s])
    dat$strat <- paste0(dat$strat, dat[, s])
  }

  phI_tab <- merge(phI_tab, des_long)

  V <- vector()
  for (i in 1:nrow(phI_tab)) {
    if (!is.null(wave1_Validated)) {
      V <- append(V, sample(x = which(dat[, "strat"] == phI_tab[i, "strat"] & dat[, wave1_Validated] == 0), size = phI_tab[i, "target"], replace = F))
    } else {
      V <- append(V, sample(x = which(dat[, "strat"] == phI_tab[i, "strat"]), size = phI_tab[i, "target"], replace = F))
    }
  }

  V <- as.numeric(seq(1, phI) %in% V[order(V)])
  return(V)
}
