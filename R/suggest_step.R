#' Recommends next audit step size for iterative grid search.
#' @name suggest_step
#' @param phII Phase II sample size.
#' @param min_n Minimum stratum size to be sampled.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 variables can be accommodated.
#' @param prev_grid_des If grid > 1, the audit from the previous iteration that was optimal.
#' @param prev_grid_delta If grid > 1, the step size from the previous iteration.
#' @param max_grid_size Integer maxium for the largest grids that will be searched. DEFAULT = 100000.
#' @return An integer.
#' @export
suggest_step <- function(phII, phI_strat, min_n, sample_on, prev_grid_des = NULL, prev_delta = NULL, max_grid_size = 500) {
  # Since each of the strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - num_strat x min_n
  num_strat <- 2 ^ length(sample_on)
  phi <- phII - num_strat * min_n

  # Suggest audit steps that are factors of each other to ensure valid grids
  steps <- vector()
  gcd <- max(seq(1, (phi - 1))[phi %% seq(1, (phi - 1)) == 0])
  while(gcd > 1) {
    steps <- append(steps, gcd)
    gcd <- max(seq(1, (gcd - 1))[gcd %% seq(1, (gcd - 1)) == 0])
  }
  steps <- c(steps, 1)
  first_step <- min(steps[sapply(X = steps, FUN = grid_size, phi = phi, num_strat = num_strat, phI_strat = phI_strat, prev_grid_des = NULL, prev_delta = NULL) < max_grid_size])

  if (is.null(prev_grid_des)) {
    return(first_step)
  } else if (length(steps) > 1) {
    steps <- c(first_step, steps[steps < first_step])
    keep <- rep(TRUE, length(steps))
    for (i in 2:length(steps)) {
      size <- grid_size(delta = steps[i], phi = phi, num_strat = num_strat, phI_strat = phI_strat, prev_grid_des = prev_grid_des, prev_delta = prev_delta)
      if (size > max_grid_size) {
        keep[i] <- FALSE
      }
    }
    return(min(steps[keep]))
  }
}
