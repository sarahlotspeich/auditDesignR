#' Adaptive grid search optimization to find the optimal design (optMLE)
#' @name optMLE_grid
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param phI_strat Phase I stratum sample sizes, named list: \code{N00, N01, N10, N11}.
#' @param min_n Minimum stratum size to be sampled.
#' @param window_mult Window multiplier for how wide the grid should look on either side of the previous iteration's optimal design. The multiplier applies to the \code{audit_steps}; for example \code{window_mult = 1} (the DEFAULT) will look at a window of 1 times each step around the previous grid's optimal design.
#' @param audit_steps A numeric vector (on the person scale). An audit step of 40 defines a grid with 40-person increments between designs. The elements of \code{audit_steps} should share common multiplies to avoid empty grid returns.
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index).
#' @param Y_val Column with the validated outcome (can be name or numeric index).
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index).
#' @param X_val Column(s) with the validated predictors (can be name or numeric index).
#' @param indiv_score Matrix of score vectors for all parameters, \code{beta} and \code{eta}.
#' @param errors_in Measurement error setting, options are \code{"Both"}, \code{"Y only"}, \code{"X only"}. Default is \code{"Both"}.
#' @param return_full_grid If \code{TRUE}, all audits from all iterations of the grid search will be return. Default is \code{FALSE}.
#' @return Scalar function value.
#' @export
optMLE_grid <- function(phI, phII, phI_strat, min_n, window_mult = 1, audit_steps = c(10, 1),
                        Y_unval, Y_val, X_unval, X_val, indiv_score, errors_in, return_full_grid = FALSE) {

  audit_windows <- window_mult * c(NA, audit_steps[-length(audit_steps)])

  # Since each of the 4 strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - 4 x min_n
  n_to_allocate <- phII - 4 * min_n

  # Calculate all_prop_inc based on step sizes all_n_inc and n_to_allocate
  audit_steps_prop <- audit_steps / n_to_allocate

  # Create a data frame to save the optimal designs from each grid search
  all_opt_des <- data.frame(grid = 1:length(audit_steps),
                            n_inc = audit_steps,
                            prop_inc = audit_steps_prop,
                            n00 = NA, n01 = NA, n10 = NA, n11 = NA,
                            pi00 = NA, pi01 = NA, pi10 = NA, pi11 = NA,
                            Vbeta = NA, grid_size = NA)

  # Run initial grid search
  new_grid <- prop_grid(prop_min = c(0, 0, 0, 0),
                        prop_max = c(1, 1, 1, 1),
                        prop_inc = audit_steps_prop[1],
                        phI_strat = phI_strat,
                        phII = phII,
                        min_n = min_n,
                        n_to_allocate = n_to_allocate)

  new_grid_list <- split(new_grid, seq(nrow(new_grid)))
  for (r in 1:length(new_grid_list)) {
    new_grid_list[[r]] <- as.numeric((new_grid_list[[r]][,c("pi00", "pi01", "pi10", "pi11")]))
  }

  new_grid$Vbeta <- sapply(X = new_grid_list,
                           FUN = var_formula,
                           Y_unval = Y_unval,
                           Y_val = Y_val,
                           X_unval = X_unval,
                           X_val = X_val,
                           phI = phI,
                           indiv_score = indiv_score,
                           errors_in = errors_in)

  min_var <- min(new_grid$Vbeta)
  # Check for a clear minimum
  if (sum(new_grid$Vbeta == min_var) == 1) {
    min_var_design <- new_grid[new_grid$Vbeta == min_var,]
    # Sort grid by V(beta) smallest -> largest
    new_grid <- new_grid[order(new_grid$Vbeta),]
    prev_min <- new_grid[1, c("n00", "n01", "n10", "n11", "pi00", "pi01", "pi10", "pi11", "Vbeta")]
    all_opt_des[1, c("n00", "n01", "n10", "n11", "pi00", "pi01", "pi10", "pi11", "Vbeta")] <-  prev_min
    all_opt_des[1, "grid_size"] <- nrow(new_grid)
  } else {
    warning("Unable to find clear minimum. Please select a new starting value.")
    return(list("all_opt" = NA,
                "min_var" = 9999,
                "min_var_design" = NA,
                "findOptimal" = FALSE,
                "full_grid_search" = NA))
  }

  if (return_full_grid) { all_grids <- cbind(grid = 1, new_grid) }

  findOptimal <- findFinalOptimal <- FALSE
  skippedLast <- FALSE
  for (step in 2:length(audit_steps_prop)) {
    # To choose min/max for next grid, use previous grid's "optimal" design
    prev_grid_allocated <- prev_min[, c("n00", "n01", "n10", "n11")] - min_n

    # If skipped the last grid, keep the wider window
    if (skippedLast) {
      # Run new grid
      new_grid <- prop_grid(prop_min = pmax(0, (prev_grid_allocated - audit_windows[step - 1]) / n_to_allocate),
                            prop_max = pmin(1, (prev_grid_allocated + audit_windows[step - 1]) / n_to_allocate),
                            prop_inc = audit_steps_prop[step],
                            phII = phII,
                            phI_strat = phI_strat,
                            min_n = min_n,
                            n_to_allocate = n_to_allocate)
    } else {
      # Run new grid
      new_grid <- prop_grid(prop_min = pmax(0, (prev_grid_allocated - audit_windows[step]) / n_to_allocate),
                            prop_max = pmin(1, (prev_grid_allocated + audit_windows[step]) / n_to_allocate),
                            prop_inc = audit_steps_prop[step],
                            phII = phII,
                            phI_strat = phI_strat,
                            min_n = min_n,
                            n_to_allocate = n_to_allocate)
    }

    if (nrow(new_grid) == 0 & step == length(audit_steps_prop)) {
      # Run new grid
      new_grid <- prop_grid(prop_min = pmax(0, (prev_grid_allocated - audit_windows[step - 1]) / n_to_allocate),
                            prop_max = pmin(1, (prev_grid_allocated + audit_windows[step - 1]) / n_to_allocate),
                            prop_inc = audit_steps_prop[step],
                            phII = phII,
                            phI_strat = phI_strat,
                            min_n = min_n,
                            n_to_allocate = n_to_allocate)
    }

    if (nrow(new_grid) > 0) {
      new_grid_list <- split(new_grid, seq(nrow(new_grid)))
      for (r in 1:length(new_grid_list)) {
        new_grid_list[[r]] <- as.numeric((new_grid_list[[r]][,c("pi00", "pi01", "pi10", "pi11")]))
      }

      new_grid$Vbeta <- sapply(X = new_grid_list,
                               FUN = var_formula,
                               Y_unval = Y_unval,
                               Y_val = Y_val,
                               X_unval = X_unval,
                               X_val = X_val,
                               phI = phI,
                               indiv_score = indiv_score,
                               errors_in = errors_in)

      # Check that a clear minimum was found
      ## And that the new minimum variance is <= the previous
      min_var <- min(new_grid$Vbeta)
      findOptimal <- sum(new_grid$Vbeta == min_var) == 1

      if (findOptimal & min_var <= all_opt_des$Vbeta[step - 1]) {
        min_var_design <- new_grid[new_grid$Vbeta == min_var, ]
        all_opt_des[step, c("n00", "n01", "n10", "n11", "pi00", "pi01", "pi10", "pi11", "Vbeta")] <- min_var_design
        all_opt_des[step, "grid_size"] <- nrow(new_grid)
        if (return_full_grid) { all_grids <- rbind(all_grids, cbind(grid = step, new_grid)) }
        if (step == length(audit_steps_prop)) { findFinalOptimal <- TRUE }
        prev_min <- min_var_design
        skippedLast <- FALSE
      } else { skippedLast <- TRUE }
    } else { skippedLast <- TRUE }
    if (skippedLast) { all_opt_des[step, "Vbeta"] <- all_opt_des[step - 1, "Vbeta"]}
  }

  if(findFinalOptimal) {
    opt_des <- all_opt_des[nrow(all_opt_des), ]
  } else {
    opt_des <- data.frame(n00 = NA, n01 = NA, n10 = NA, n11 = NA)
  }
  if (!return_full_grid) { all_grids = NA }
  return(list("all_opt" = all_opt_des,
              "min_var" = min_var,
              "min_var_design" = min_var_design,
              "findOptimal" = findFinalOptimal,
              "full_grid_search" = all_grids))
}
