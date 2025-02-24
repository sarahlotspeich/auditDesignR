#' Sample unvalidated case-control (CC*) audit based on Phase I variables.
#' @name sample_cc
#' @param dat Dataframe or matrix containing columns \code{sample_on}.
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param sample_on Column with the Phase I outcome variable (should be categorical) used for sampling strata (can be name or numeric index).
#' @param prop_cases Proportion of the audit to be taken from cases (based on \code{sample_on}). Default is \code{prop_cases = 0.5} for even numbers of cases and controls.
#' @param wave1_Validated (For use with multi-wave designs) A logical vector with the validation indicator from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_cc <- function(dat, phI, phII, sample_on, prop_cases = 0.5, wave1_Validated = NULL) {
  if (length(wave1_Validated) == 0) {
    wave1_Validated = rep(FALSE, phI)
  }
  
  # Cross-tabulate Phase I data (Y*)
  phI_tab <- data.frame(table(dat[!wave1_Validated, sample_on]))
  colnames(phI_tab)[1] <- sample_on
  colnames(phI_tab)[2] <- "N"
  phI_tab$target <- pmin(floor(phII * c(1 - prop_cases, prop_cases)), phI_tab$N)

  # If any of the strata were too small, redistribute leftover
  while(sum(phI_tab$target) < phII) {
    for (i in 1:nrow(phI_tab)) {
      if (phI_tab$N[i] > phI_tab$target[i] & sum(phI_tab$target) < phII) {
        phI_tab$target[i] <- phI_tab$target[i] + 1
      }
    }
  }

  V_cc <- vector()
  for (i in 1:nrow(phI_tab)) {
    V_cc <- append(V_cc, 
                   sample(x = which(dat[, sample_on] == phI_tab[i, sample_on] & !wave1_Validated), 
                          size = phI_tab[i, "target"], 
                          replace = F))
  }

  V_cc <- as.numeric(seq(1, phI) %in% V_cc[order(V_cc)])
  return(V_cc)
}
