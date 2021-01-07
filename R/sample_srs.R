#' Sample simple random sample (SRS) audit.
#' @name sample_srs
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @return A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
#' @export
sample_srs <- function(phI, phII) {
  V_srs <- as.numeric(seq(1, phI) %in% sample(x = seq(1, phI), size = phII, replace = F))
  return(V_srs)
}
