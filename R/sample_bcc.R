#' Sample unvalidated balanced case-control (BCC*) audit based on Phase I variables.
#' @name sample_bcc
#' @param dat Dataframe or matrix containing columns \code{sample_on}.
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 binary variables can be accommodated.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_bcc <- function(dat, phI, phII, sample_on) {

  # Cross-tabulate Phase I data (Y*)
  phI_tab <- data.frame(table(dat[, sample_on]))
  colnames(phI_tab)[1:length(sample_on)] <- sample_on
  colnames(phI_tab)[ncol(phI_tab)] <- "N"
  phI_tab$target <- pmin(floor(phII / nrow(phI_tab)), phI_tab$N)

  # Create stratum IDs to sample on
  phI_tab$strat <- "N"
  dat$strat <- "N"
  for (s in sample_on) {
    phI_tab$strat <- paste0(phI_tab$strat, phI_tab[, s])
    dat$strat <- paste0(dat$strat, dat[, s])
  }

  # If any of the strata were too small, redistribute leftover
  while(sum(phI_tab$target) < phII) {
    for (i in 1:nrow(phI_tab)) {
      if (phI_tab$N[i] > phI_tab$target[i] & sum(phI_tab$target) < phII) {
        phI_tab$target[i] <- phI_tab$target[i] + 1
      }
    }
  }

  V <- vector()
  for (i in 1:nrow(phI_tab)) {
    V <- append(V, sample(x = which(dat[, "strat"] == phI_tab[i, "strat"]), size = phI_tab[i, "target"], replace = F))
  }

  V <- as.numeric(seq(1, phI) %in% V[order(V)])
  return(V)
}
