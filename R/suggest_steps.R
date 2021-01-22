#' Adaptive grid search optimization to find the optimal design (optMLE)
#' @name suggest_steps
#' @param phII Phase II sample size.
#' @param min_n Minimum stratum size to be sampled.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 variables can be accommodated.
#' @return A vector of integers.
#' @export
suggest_steps <- function(phII, min_n, sample_on) {
  # Since each of the strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - num_strat x min_n
  num_strat <- 2 ^ length(sample_on)
  n_to_allocate <- phII - num_strat * min_n

  # Suggest audit steps that are factors of each other to ensure valid grids
  steps <- vector()
  gcd <- max(seq(1, (n_to_allocate - 1))[n_to_allocate %% seq(1, (n_to_allocate - 1)) == 0])
  while(gcd > 1) {
    steps <- append(steps, gcd)
    gcd <- max(seq(1, (gcd - 1))[gcd %% seq(1, (gcd - 1)) == 0])
  }

  return(c(steps, 1))
}
