#' Adaptive grid search optimization to find the optimal design (optMLE)
#' @name optMLE_grid
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param phI_strat Phase I stratum sample sizes, named list: \code{N00, N01, N10, N11}.
#' @param min_n Minimum stratum size to be sampled.
#' @param window_mult Window multiplier for how wide the grid should look on either side of the previous iteration's optimal design. The multiplier applies to the \code{audit_steps}; for example \code{window_mult = 1} (the DEFAULT) will look at a window of 1 times each step around the previous grid's optimal design.
#' @param audit_steps A numeric vector (on the person scale). An audit step of 40 defines a grid with 40-person increments between designs. The elements of \code{audit_steps} should share common multiplies to avoid empty grid returns.
#' @param max_window_mult If \code{window_mult} fails to find a minimum, it will increase up to \code{max_window_mult} and try again. The DEFAULT is \code{max_window_mult = 3}.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 variables can be accommodated.
#' @param indiv_score Matrix of score vectors for all parameters, \code{beta} and \code{eta}.
#' @param return_full_grid If \code{TRUE}, all audits from all iterations of the grid search will be return. Default is \code{FALSE}.
#' @return Scalar function value.
#' @export
optMLE_grid <- function(phI, phII, phI_strat, min_n, window_mult = 1, audit_steps, max_window_mult = 3,
                        sample_on, indiv_score, return_full_grid = FALSE) {
  #audit_windows <- window_mult * c(NA, audit_steps[-length(audit_steps)])
  audit_windows <- ceiling(c(NA, audit_steps[-length(audit_steps)]) / audit_steps ) * audit_steps #c(NA, audit_steps[-length(audit_steps)])
  num_strat <- 2 ^ length(sample_on)
  # Since each of the strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - num_strat x min_n
  n_to_allocate <- phII - num_strat * min_n
  audit_steps_prop <- audit_steps / n_to_allocate
  if (length(sample_on) == 1) {
    # Create a data frame to save the optimal designs from each grid search
    all_opt_des <- data.frame(grid = 1:length(audit_steps),
                              audit_steps, audit_windows,
                              prop_inc = audit_steps_prop,
                              n0 = NA, n1 = NA, pi0 = NA, pi1 = NA,
                              Vbeta = NA, grid_size = NA)
  } else if (length(sample_on) == 2) {
    # Create a data frame to save the optimal designs from each grid search
    all_opt_des <- data.frame(grid = 1:length(audit_steps),
                              audit_steps, audit_windows,
                              prop_inc = audit_steps_prop,
                              n00 = NA, n01 = NA, n10 = NA, n11 = NA,
                              pi00 = NA, pi01 = NA, pi10 = NA, pi11 = NA,
                              Vbeta = NA, grid_size = NA)
  } else if (length(sample_on) == 3) {
    # Create a data frame to save the optimal designs from each grid search
    all_opt_des <- data.frame(grid = 1:length(audit_steps),
                              audit_steps, audit_windows,
                              prop_inc = audit_steps_prop,
                              n00_0 = NA, n01_0 = NA, n10_0 = NA, n11_0 = NA,
                              n00_1 = NA, n01_1 = NA, n10_1 = NA, n11_1 = NA,
                              pi00_0 = NA, pi01_0 = NA, pi10_0 = NA, pi11_0 = NA,
                              pi00_1 = NA, pi01_1 = NA, pi10_1 = NA, pi11_1 = NA,
                              Vbeta = NA, grid_size = NA)
  } else {
    return(warning("You are attempting to sample on too many variables."))
  }

  # Run initial grid search
  new_grid <- prop_grid(prop_min = rep(0, num_strat),
                        prop_max = rep(1, num_strat),
                        prop_inc = audit_steps_prop[1],
                        phI_strat = phI_strat,
                        phII = phII,
                        min_n = min_n,
                        n_to_allocate = n_to_allocate)

  if (nrow(new_grid) == 0) {
    warning("Insufficient grid space. Please select a new starting value.")
    return(list("all_opt" = NA,
                "min_var" = 9999,
                "min_var_design" = NA,
                "findOptimal" = FALSE,
                "full_grid_search" = NA,
                "message" = "Insufficient starting grid size"))
  }

  new_grid_list <- split(new_grid, seq(nrow(new_grid)))
  for (r in 1:length(new_grid_list)) {
    new_grid_list[[r]] <- as.numeric((new_grid_list[[r]][, grep(pattern = "pi", colnames(new_grid), value = TRUE)]))
  }

  new_grid$Vbeta <- sapply(X = new_grid_list,
                           FUN = var_formula,
                           phI = phI,
                           indiv_score = indiv_score,
                           sample_on = sample_on)

  min_var <- min(new_grid$Vbeta)

  if (is.na(min_var)) {
    return(list("all_opt" = NA,
                "min_var" = 9999,
                "min_var_design" = NA,
                "findOptimal" = FALSE,
                "full_grid_search" = NA,
                "message" = "Singular information"))
  }

  # Check for a clear minimum
  findOptimal <- sum(new_grid$Vbeta == min_var) == 1
  if (findOptimal) {
    min_var_design <- prev_min <- new_grid[new_grid$Vbeta == min_var,]
    all_opt_des[1, -c(1:4)] <-  c(min_var_design, nrow(new_grid))
  } else {
    warning("Unable to find clear minimum. Please select a new starting value.")
    return(list("all_opt" = NA,
                "min_var" = 9999,
                "min_var_design" = NA,
                "findOptimal" = FALSE,
                "full_grid_search" = NA,
                "message" = "Insufficient starting grid size"))
  }

  if (return_full_grid) { all_grids <- cbind(grid = 1, new_grid) }

  #findOptimal <- findFinalOptimal <- FALSE
  for (step in 2:length(audit_steps_prop)) {
    # Reset indicator
    findOptimal <- findFinalOptimal <- FALSE

    # Reset audit window for this step
    step_mult <- window_mult

    # To choose min/max for next grid, use previous grid's "optimal" design
    prev_grid_allocated <- prev_min[, 1:num_strat] - min_n

    # Run new grid
    new_grid <- prop_grid(prop_min = pmax(0, (prev_grid_allocated - audit_windows[step]) / n_to_allocate),
                          prop_max = pmin(1, (prev_grid_allocated + audit_windows[step]) / n_to_allocate),
                          prop_inc = audit_steps_prop[step],
                          prev_grid = prev_grid_allocated,
                          phII = phII,
                          phI_strat = phI_strat,
                          min_n = min_n,
                          n_to_allocate = n_to_allocate)

    # If nrow(new_grid) = 0, widen window
    while (nrow(new_grid) == 0 & step_mult < max_window_mult) {
      step_mult <- step_mult + 1
      audit_windows[step] <- step_mult * c(NA, audit_steps[-length(audit_steps)])[step]
      new_grid <- prop_grid(prop_min = pmax(0, (prev_grid_allocated - audit_windows[step]) / n_to_allocate),
                            prop_max = pmin(1, (prev_grid_allocated + audit_windows[step]) / n_to_allocate),
                            prop_inc = audit_steps_prop[step],
                            phII = phII,
                            phI_strat = phI_strat,
                            min_n = min_n,
                            n_to_allocate = n_to_allocate)
    }

    if (nrow(new_grid) > 0) {
      new_grid_list <- split(new_grid, seq(nrow(new_grid)))
      for (r in 1:length(new_grid_list)) {
        new_grid_list[[r]] <- as.numeric((new_grid_list[[r]][, grep(pattern = "pi", colnames(new_grid), value = TRUE)]))
      }

      new_grid$Vbeta <- sapply(X = new_grid_list,
                               FUN = var_formula,
                               phI = phI,
                               indiv_score = indiv_score,
                               sample_on = sample_on)

      # Check that a clear minimum was found
      ## And that the new minimum variance is <= the previous
      min_var <- min(new_grid$Vbeta)

      if (is.na(min_var)) {
        return(list("all_opt" = NA,
                    "min_var" = 9999,
                    "min_var_design" = NA,
                    "findOptimal" = FALSE,
                    "full_grid_search" = NA,
                    "message" = "Singular information"))
      }

      findOptimal <- sum(new_grid$Vbeta == min_var) == 1 & (min_var - all_opt_des$Vbeta[step - 1]) < 1E-8

      if (findOptimal) {
        min_var_design <- prev_min <- new_grid[new_grid$Vbeta == min_var,]
        all_opt_des[step, -c(1:4)] <-  c(min_var_design, nrow(new_grid))
        if (return_full_grid) { all_grids <- rbind(all_grids, cbind(grid = step, new_grid)) }
        if (step == length(audit_steps_prop)) { findFinalOptimal <- TRUE }
        try_wider <- FALSE
      } else { try_wider <- TRUE }
    }

    if (try_wider) {
      while ((min_var - all_opt_des$Vbeta[step - 1]) >= 1E-8 & step_mult < max_window_mult) {
        step_mult <- step_mult + 1
        audit_windows[step] <- step_mult * c(NA, audit_steps[-length(audit_steps)])[step]
        new_grid <- prop_grid(prop_min = pmax(0, (prev_grid_allocated - audit_windows[step]) / n_to_allocate),
                              prop_max = pmin(1, (prev_grid_allocated + audit_windows[step]) / n_to_allocate),
                              prop_inc = audit_steps_prop[step],
                              phII = phII,
                              phI_strat = phI_strat,
                              min_n = min_n,
                              n_to_allocate = n_to_allocate)

        if (nrow(new_grid) > 0) {
          new_grid_list <- split(new_grid, seq(nrow(new_grid)))
          for (r in 1:length(new_grid_list)) {
            new_grid_list[[r]] <- as.numeric((new_grid_list[[r]][, grep(pattern = "pi", colnames(new_grid), value = TRUE)]))
          }

          new_grid$Vbeta <- sapply(X = new_grid_list,
                                   FUN = var_formula,
                                   phI = phI,
                                   indiv_score = indiv_score,
                                   sample_on = sample_on)

          # Check that a clear minimum was found
          ## And that the new minimum variance is <= the previous
          min_var <- min(new_grid$Vbeta)

          if (is.na(min_var)) {
            return(list("all_opt" = NA,
                        "min_var" = 9999,
                        "min_var_design" = NA,
                        "findOptimal" = FALSE,
                        "full_grid_search" = NA,
                        "message" = "Singular information"))
          }

          findOptimal <- sum(new_grid$Vbeta == min_var) == 1 & (min_var - all_opt_des$Vbeta[step - 1]) < 1E-8

          if (findOptimal) {
            min_var_design <- prev_min <- new_grid[new_grid$Vbeta == min_var,]
            all_opt_des[step, -c(1:4)] <-  c(min_var_design, nrow(new_grid))
            if (return_full_grid) { all_grids <- rbind(all_grids, cbind(grid = step, new_grid)) }
            if (step == length(audit_steps_prop)) { findFinalOptimal <- TRUE }
          }
        }
      }
    }

    if (!findOptimal) {
      all_opt_des[step, "Vbeta"] <- all_opt_des[step - 1, "Vbeta"]
    }
  }

  all_opt_des$audit_windows <- audit_windows

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
              "full_grid_search" = all_grids,
              "message" = ifelse(findFinalOptimal, "Grid search successful", "Grid completed without finding minimum")))
}
