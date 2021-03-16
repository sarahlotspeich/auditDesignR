#' Recommends next audit step size for iterative grid search.
#' @name suggest_step
#' @param phII Phase II sample size.
#' @param min_n Minimum stratum size to be sampled.
#' @param num_strat Number of strata on which sampling is based.
#' @param num_steps Character option for how many steps to be returned, options include \code{"smallest"} and \code{"all"}. DEFAULT is \code{"smallest"}.
#' @param prev_grid_des If grid > 1, the audit from the previous iteration that was optimal.
#' @param prev_grid_delta If grid > 1, the step size from the previous iteration.
#' @param max_grid_size Integer maxium for the largest grids that will be searched.
#' @return An integer.
#' @export
suggest_step <- function(phII, phI_strat, min_n, num_strat, num_steps = "smallest", prev_grid_des, prev_delta, max_grid_size) {
  # Since each of the strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - num_strat x min_n
  phi <- phII - sum(unlist(pmin(phI_strat, min_n)))

  # Suggest audit steps that are factors of each other to ensure valid grids
  steps <- vector()
  gcd <- max(seq(1, (phi - 1))[phi %% seq(1, (phi - 1)) == 0])
  while(gcd > 1) {
    steps <- append(steps, gcd)
    gcd <- max(seq(1, (gcd - 1))[gcd %% seq(1, (gcd - 1)) == 0])
  }
  steps <- c(steps, 1)
  small_enough <- steps[sapply(X = steps, FUN = grid_size, phi = phi, num_strat = num_strat, phI_strat = phI_strat, prev_grid_des = NULL, prev_delta = NULL) < max_grid_size]
  first_step <- min(small_enough)

  if (is.null(prev_grid_des)) {
    if (num_steps == "smallest") {
      return(first_step)
    } else if (num_steps == "all") {
      return(steps)
    }
  } else if (length(steps) > 1 & num_steps == "smallest") {
    if (first_step == prev_delta) {
      steps <- c(steps[steps < prev_delta])
    } else {
      steps <- c(prev_delta, steps[steps < prev_delta])
    }
    steps <- c(first_step, steps)
    keep <- rep(TRUE, length(steps))
    keep[steps >= prev_delta] <- FALSE
    for (i in 2:length(steps)) {
      size <- grid_size(delta = steps[i], phi = phi, num_strat = num_strat, phI_strat = phI_strat, prev_grid_des = prev_grid_des, prev_delta = prev_delta)
      if (size > max_grid_size | size <= 1) {
        keep[i] <- FALSE
      }
    }
    # If no factors will work, consider all numbers smaller (if previous was small enough)
    if (mean(keep) == 0 & prev_delta <= 20) {
      smaller_steps <- seq(1, (prev_delta - 1))
      keep <- rep(TRUE, length(smaller_steps))
      for (i in 1:length(smaller_steps)) {
        size <- grid_size(delta = smaller_steps[i], phi = phi, num_strat = num_strat, phI_strat = phI_strat, prev_grid_des = prev_grid_des, prev_delta = prev_delta)
        if (size > max_grid_size) {
          keep[i] <- FALSE
        }
      }
      if (mean(keep) == 0) {
        return(9999)
      } else {
        return(min(smaller_steps[keep]))
      }
    }

    if (mean(keep) == 0) {
      return(9999)
    } else {
      return(min(steps[keep]))
    }
  }
}
