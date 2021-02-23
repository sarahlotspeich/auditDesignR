#' Build audit grid
#' @name grid_size
#' @param delta Audit step size (in people).
#' @param phi Number of people to be allocated.
#' @param num_strat Number of strata on which sampling is based. Currently handles \code{num_strat} = 2, 4, or 8.
#' @param phI_strat Phase I stratum sample sizes as a named list, dataframe, or vector
#' @param phIIa_strat For multi-wave designs, Phase II(a) stratum sample sizes as a named list, dataframe, or vector. Default is \code{NULL}.
#' @param closed For multi-wave designs, a vector of names for strata that are "closed", meaning we do not wish to sample from them. Default is \code{NULL}.
#' @param closed_at For multi-wave designs, a vector of already sampled sizes for strata that are "closed" (must be the same length as \code{closed}). Default is \code{NULL}.
#' @param min_n Minimum stratum size to be sampled.
#' @param prev_grid_des If grid > 1, the audit from the previous iteration that was optimal.
#' @param prev_delta If grid > 1, the step size from the previous iteration.
#' @export
build_grid <- function(delta, phi, num_strat, phI_strat, phIIa_strat = NULL, closed = NULL, closed_at = NULL, min_n, prev_grid_des = NULL, prev_delta = NULL) {
  stars <- round(phi / delta)

  # Exclude closed strata
  if (length(closed) > 0) {
    phI_strat <- phI_strat[names(phI_strat) !=  toupper(closed)]
  }

  # If Phase II(a) strata were provided, these need to be substracted from phI_strat
  ## Since they are no longer available for audit in Phase II(b)
  if (!is.null(phIIa_strat)) {
    phI_strat_rem <- phI_strat
    for (s in names(phI_strat_rem)) {
      phI_strat_rem[s] <- phI_strat_rem[s] - phIIa_strat[tolower(s)]
    }
  }

  # First grid tries entire space
  if (is.null(prev_grid_des)) {
    window_lb <- rep(0, num_strat)
    if (!is.null(phIIa_strat)) {
      window_ub <- pmin(stars, floor((unlist(phI_strat_rem) - pmin(unlist(phI_strat_rem), min_n)) / delta))
    } else {
      window_ub <- pmin(stars, floor((unlist(phI_strat) - pmin(unlist(phI_strat), min_n)) / delta))
    }
  } else {
    if (!is.null(phIIa_strat)) {
      window_lb <- floor(pmax(prev_grid_des - prev_delta, rep(0, num_strat)) / delta)
      window_ub <- floor(pmin(prev_grid_des + prev_delta, (unlist(phI_strat_rem) - pmin(unlist(phI_strat_rem), min_n))) / delta)
    } else {
      window_lb <- floor(pmax(prev_grid_des - prev_delta, rep(0, num_strat)) / delta)
      window_ub <- floor(pmin(prev_grid_des + prev_delta, (unlist(phI_strat) - pmin(unlist(phI_strat), min_n))) / delta)
    }
  }

  grid <- do.call(expand.grid, grid_vals(x = window_lb, y = window_ub))

  # include previous if not included
  if (!is.null(prev_grid_des)) {
    grid <- unique(rbind(grid, prev_grid_des / delta))
  }

  # filter to grids that sum to audit constraint
  grid <- grid[rowSums(grid) == stars, ]

  # scale back to people
  grid <- grid * delta

  # add back min_n to each stratum
  grid <- grid +
    matrix(data = pmin(unlist(phI_strat), min_n), nrow = nrow(grid), ncol = ncol(grid), byrow = TRUE)
  colnames(grid) <- gsub("N", "n", names(phI_strat))

  # If multi-wave design, also add the stratum sample sizes from Phase II(a)
  if (!is.null(phIIa_strat)) {
    for (c in colnames(grid)) {
      grid[, c] <- grid[, c] + as.numeric(phIIa_strat[c])
    }
  }

  # Check for grids where n > N (truncate at N)
  for (c in 1:ncol(grid)) {
    grid[, c]
  }

  # Augment with sampling probabilities
  pi_grid <- matrix(nrow = nrow(grid), ncol = ncol(grid))
  colnames(pi_grid) <- gsub("n", "pi", colnames(grid))

  for (c in 1:ncol(pi_grid)) {
    pi_grid[, c] <- grid[, c] / phI_strat[[c]]
  }

  grid <- cbind(grid, pi_grid)

  return(grid)
}
