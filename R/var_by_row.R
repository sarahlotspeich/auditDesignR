#' @export
var_by_row <- function(row, phI, indiv_score, sample_on) {
  pi_vec <- row[(length(row)/2+1):length(row)]
  Vbeta <- var_formula(pi_vec = pi_vec, phI = phI, indiv_score = indiv_score, sample_on = sample_on)
  return(Vbeta)
}
