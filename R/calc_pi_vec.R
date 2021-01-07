#' Calculate the validation/sampling probabilities for a two-phase design.
#' @name calc_pi_vec
#' @param Validated Columns with the validation indicator (can be name or numeric index)
#' @return 1 x 4 vector of the sampling probabilities, pi:
#' \item{pi00}{P(V=1|\code{phI_Y}=0, \code{phI_X}=0)}
#' \item{pi01}{P(V=1|\code{phI_Y}=0, \code{phI_X}=1)}
#' \item{pi10}{P(V=1|\code{phI_Y}=1, \code{phI_X}=0)}
#' \item{pi11}{P(V=1|\code{phI_Y}=1, \code{phI_X}=1)}
#' @export
calc_pi_vec <- function(Validated, dat, phI_Y, phI_X) {
  dat_val <- dat[Validated == 1, ]
  tab_phI <- table(dat[, phI_Y], dat[, phI_X])
  tab_phII <- table(dat_val[, phI_Y], dat_val[, phI_X])
  pi_vec = c(pi00 = tab_phII[1, 1] / tab_phI[1, 1], pi01 = tab_phII[1, 2] / tab_phI[1, 2],
             pi10 = tab_phII[2, 1] / tab_phI[2, 1], pi11 = tab_phII[2, 2] / tab_phI[2, 2])
  return(pi_vec)
}
