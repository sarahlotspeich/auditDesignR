#' Sample simple random sample (SRS) audit.
#' @name sample_srs
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param wave1_Validated (For use with multi-wave designs) A logical vector with the validation indicator from first wave of audits.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_srs <- function(phI, phII, wave1_Validated = NULL) {
  if (length(wave1_Validated) == 0) {
    wave1_Validated = rep(FALSE, phI)
  }
  V_srs <- as.numeric(seq(1, phI) %in% sample(x = seq(1, phI)[!wave1_Validated], size = phII, replace = F))
  return(V_srs)
}
