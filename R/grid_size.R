#' Calculate size of the audit grid
#' @name grid_size
#' @param delta Audit step size (in people).
#' @param phi Number of people to be allocated: \code{phi = phII - num_strat * min_n}.
#' @param num_strat Number of strata on which sampling is based. Currently handles \code{num_strat} = 2, 4, or 8.
#' @param phI_strat Phase I stratum sample sizes, named list.
#' @param prev_grid_des If grid > 1, the audit from the previous iteration that was optimal.
#' @param prev_delta If grid > 1, the step size from the previous iteration.
#' @return An integer.
grid_size <- function(delta, phi, num_strat, phI_strat, prev_grid_des = NULL, prev_delta = NULL) {
  # Stars and bars
  stars <- phi / delta
  bars <- num_strat - 1

  if (is.null(prev_grid_des)) {
    return(choose(n = (stars + bars), k = (bars)))
  } else {
    # Number of all solutions (before taking prop_min/prop_max into account)
    window_lb <- pmax(prev_grid_des - prev_delta, rep(0, num_strat)) / delta
    window_ub <- pmin(prev_grid_des + prev_delta, unlist(phI_strat)) / delta

    # Number of solutions >= window_lb
    grid_size_lb <- choose(n = (stars + bars - sum(window_lb)), k = bars)
    num_too_big <- 0

    for (s in 1:length(window_ub)) {
      too_big <- as.numeric(stars - (window_ub[s] + 1))
      while (too_big >= 0 & choose(n = (too_big + (bars - 1) - sum(window_lb[-s])), k = (bars - 1)) > 0) {
        num_too_big <- num_too_big + choose(n = (too_big + (bars - 1) - sum(window_lb[-s])), k = (bars - 1))
        too_big <- too_big - 1
      }
    }
    return(grid_size_lb - num_too_big)
  }
}
