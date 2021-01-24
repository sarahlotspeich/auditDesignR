#' Build audit grid
#' @name grid_size
#' @param delta Audit step size (in people).
#' @param phi Number of people to be allocated: \code{phi = phII - num_strat * min_n}.
#' @param num_strat Number of strata on which sampling is based. Currently handles \code{num_strat} = 2, 4, or 8.
#' @param phI_strat Phase I stratum sample sizes, named list.
#' @param min_n Minimum stratum size to be sampled.
#' @param prev_grid_des If grid > 1, the audit from the previous iteration that was optimal.
#' @param prev_delta If grid > 1, the step size from the previous iteration.
#' @export
build_grid <- function(delta, phi, num_strat, phI_strat, min_n, prev_grid_des = NULL, prev_delta = NULL) {
  stars <- phi / delta
  # First grid tries entire space
  if (is.null(prev_grid_des)) {
    if (grid_size(delta = delta, phi = phi, num_strat = num_strat, phI_strat = phI_strat, prev_grid_des = NULL, prev_delta = NULL) < 1000) {
      window_lb <- rep(0, num_strat)
      window_ub <- pmin(stars, floor(unlist(phI_strat) / delta))
    } else {
      return(warning("Initial grid is too large. Please increase delta and try again."))
    }
  } else {
    window_lb <- pmax(prev_grid_des - prev_delta, rep(0, num_strat)) / delta
    window_ub <- pmin(prev_grid_des + prev_delta, unlist(phI_strat)) / delta

  }
  grid <- do.call(expand.grid, mapply(":", window_lb, window_ub))
  grid <- grid[rowSums(grid) == stars, ]
  grid <- grid * delta + min_n
  colnames(grid) <- gsub("N", "n", names(phI_strat))

  # Augment with sampling probabilities
  pi_grid <- matrix(nrow = nrow(grid), ncol = ncol(grid))
  colnames(pi_grid) <- gsub("n", "pi", colnames(grid))

  for (c in 1:ncol(pi_grid)) {
    pi_grid[, c] <- grid[, c] / phI_strat[[c]]
  }

  grid <- cbind(grid, pi_grid)

  return(grid)
}
