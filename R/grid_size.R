#' Calculate size of the audit grid
#' @name grid_size
#' @param delta Audit step size (in people).
#' @param phi Number of people to be allocated (after \code{min_n}).
#' @param num_strat Number of strata on which sampling is based. Currently handles \code{num_strat} = 2, 4, or 8.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 variables can be accommodated.
#' @param prev_grid_des If grid > 1, the audit from the previous iteration that was optimal.
#' @return A vector of integers.
grid_size <- function(delta, phi, num_strat, prev_grid_des = NULL) {
  # Stars and bars
  stars <- phi / delta
  bars <- num_strat - 1

  if (is.null(prev_grid_des)) {
    return(choose(n = (stars + bars), k = (bars)))
  } else {
    # Number of all solutions (before taking prop_min/prop_max into account)
    # all_sol <- choose(n = (stars + bars), k = (bars))
    # Number of all solutions (before taking prop_min into account)

    # Subtract off "bad" solutions, given the prop_min/prop_max window constraint
  }
}
