run_grid <- function(new_grid, Y_unval = NULL, Y_val, X_unval = NULL, X_val, phI, indiv_score, errors_in) {
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
    status <- TRUE
  } else {
    min_var_design <- data.frame(n00 = NA, n01 = NA, n10 = NA, n11 = NA, pi00 = NA, pi01 = NA, pi10 = NA, pi11 = NA, Vbeta = NA)
    status <- FALSE
  }

  new_grid <- new_grid[order(new_grid$Vbeta),]

  return(list("grid" = new_grid, "min_var" = min_var, "min_var_design" = min_var_design, "find_optimal" = status))
}

prop_grid <- function(prop_min, prop_max, prop_inc, phI_strat, phII, min_n, n_to_allocate) {
  grid <- expand.grid(prop_n00 = seq(prop_min[1], prop_max[1], by = prop_inc),
                      prop_n01 = seq(prop_min[2], prop_max[2], by = prop_inc),
                      prop_n10 = seq(prop_min[3], prop_max[3], by = prop_inc),
                      prop_n11 = seq(prop_min[4], prop_max[4], by = prop_inc))
  grid <- grid[rowSums(grid) == 1, ]

  grid$n00 <- grid$prop_n00 * n_to_allocate #phI_strat$N00
  grid$n01 <- grid$prop_n01 * n_to_allocate #phI_strat$N01
  grid$n10 <- grid$prop_n10 * n_to_allocate #phI_strat$N10
  grid$n11 <- grid$prop_n11 * n_to_allocate #phI_strat$N11

  grid[, c("n00", "n01", "n10", "n11")] <- grid[, c("n00", "n01", "n10", "n11")] + min_n

  # Make sure n = phII
  grid[, "n"] <- round(grid[, "n00"] + grid[, "n01"] + grid[, "n10"] + grid[, "n11"])
  grid <- grid[grid[, "n"] == phII, ]

  # Make sure n < N in each stratum
  grid <- grid[grid[, "n00"] <= phI_strat$N00, ]
  grid <- grid[grid[, "n01"] <= phI_strat$N01, ]
  grid <- grid[grid[, "n10"] <= phI_strat$N10, ]
  grid <- grid[grid[, "n11"] <= phI_strat$N11, ]

  # Create columns with the P(Vi=1|Yi*, Xi*)
  grid[, "pi00"] <- grid[, "n00"] / phI_strat$N00
  grid[, "pi01"] <- grid[, "n01"] / phI_strat$N01
  grid[, "pi10"] <- grid[, "n10"] / phI_strat$N10
  grid[, "pi11"] <- grid[, "n11"] / phI_strat$N11

  return(grid[, c("n00", "n01", "n10", "n11", "pi00", "pi01", "pi10", "pi11")])
}

#' @export
optimal_audit_grid <- function(all_n_inc = c(40, 20, 10, 5, 1), window_mult = 3, phI_strat, phI, phII, min_n, indiv_score, Y_unval = NULL, Y_val, X_unval = NULL, X_val, errors_in) {
  # Since each of the 4 strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - 4 x min_n
  n_to_allocate <- phII - 4 * min_n

  # Calculate all_prop_inc based on step sizes all_n_inc and n_to_allocate
  all_prop_inc <- all_n_inc / n_to_allocate

  # Create a data frame to save the optimal designs from each grid search
  all_opt_des <- data.frame(grid = 1:length(all_n_inc),
                            n_inc = all_n_inc,
                            prop_inc = all_prop_inc,
                            n00 = NA, n01 = NA, n10 = NA, n11 = NA,
                            pi00 = NA, pi01 = NA, pi10 = NA, pi11 = NA,
                            Vbeta = NA)

  # Run initial grid search
  prev_grid <- prop_grid(prop_min = c(0, 0, 0, 0),
                         prop_max = c(1, 1, 1, 1),
                         prop_inc = all_prop_inc[1],
                         phI_strat = phI_strat,
                         phII = phII,
                         min_n = min_n,
                         n_to_allocate = n_to_allocate)

  prev_grid %>% run_grid(Y_unval = Y_unval, Y_val = Y_val, X_unval = X_unval, X_val = X_val, phI = phI, indiv_score = indiv_score, errors_in = errors_in) -> prev_grid_res
  all_opt_des[1, c("n00", "n01", "n10", "n11", "pi00", "pi01", "pi10", "pi11", "Vbeta")] <- prev_grid_res$min_var_design

  for (step in 2:length(all_prop_inc)) {
    # To choose min/max for next grid, use previous grid's "optimal" design
    prev_grid_allocated <- prev_grid_res$min_var_design[, c("n00", "n01", "n10", "n11")] - min_n
    # Run new grid
    prev_grid <- prop_grid(prop_min = pmax(0, prev_grid_allocated / n_to_allocate - window_mult * all_prop_inc[step]),
                           prop_max = pmin(1, prev_grid_allocated / n_to_allocate + window_mult * all_prop_inc[step]),
                           prop_inc = all_prop_inc[step],
                           phII = phII,
                           phI_strat = phI_strat,
                           min_n = min_n,
                           n_to_allocate = n_to_allocate)

    if (nrow(prev_grid) > 0) {
      prev_grid %>% run_grid(Y_unval = Y_unval, Y_val = Y_val, X_unval = X_unval, X_val = X_val, phI = phI, indiv_score = indiv_score, errors_in = errors_in) -> prev_grid_res
      if (prev_grid_res$find_optimal) {
        all_opt_des[step, c("n00", "n01", "n10", "n11", "pi00", "pi01", "pi10", "pi11", "Vbeta")] <- prev_grid_res$min_var_design
      }
    }
  }

  return(list("all_opt" = all_opt_des, "min_var" = prev_grid_res$min_var, "min_var_design" = prev_grid_res$min_var_design, "find_optimal" = prev_grid_res$find_optimal))
}
