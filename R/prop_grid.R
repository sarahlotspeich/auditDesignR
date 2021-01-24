prop_grid <- function(prop_min, prop_max, prop_inc, prev_grid = NULL, phI_strat, phII, min_n, n_to_allocate) {
  num_strat <- length(phI_strat)
  if (num_strat == 2) {
    grid <- expand.grid(prop_n0 = unique(seq(prop_min[1], prop_max[1], by = prop_inc), prev_grid[[1]] / n_to_allocate),
                        prop_n1 = unique(seq(prop_min[2], prop_max[2], by = prop_inc), prev_grid[[2]] / n_to_allocate))

    grid <- grid[rowSums(grid) == 1, ]

    grid$n0 <- grid$prop_n0 * n_to_allocate
    grid$n1 <- grid$prop_n1 * n_to_allocate
  } else if (num_strat == 4) {
    grid <- expand.grid(prop_n00 = unique(seq(prop_min[1], prop_max[1], by = prop_inc), prev_grid[[1]] / n_to_allocate),
                        prop_n01 = unique(seq(prop_min[2], prop_max[2], by = prop_inc), prev_grid[[2]] / n_to_allocate),
                        prop_n10 = unique(seq(prop_min[3], prop_max[3], by = prop_inc), prev_grid[[3]] / n_to_allocate),
                        prop_n11 = unique(seq(prop_min[4], prop_max[4], by = prop_inc), prev_grid[[4]] / n_to_allocate))

    grid <- grid[rowSums(grid) == 1, ]

    grid$n00 <- grid$prop_n00 * n_to_allocate
    grid$n01 <- grid$prop_n01 * n_to_allocate
    grid$n10 <- grid$prop_n10 * n_to_allocate
    grid$n11 <- grid$prop_n11 * n_to_allocate

    #grid[, c("n00", "n01", "n10", "n11")] <- grid[, c("n00", "n01", "n10", "n11")] + min_n
  } else if (num_strat == 8) {
    grid <- expand.grid(prop_n00_0 = unique(seq(prop_min[1], prop_max[1], by = prop_inc), prev_grid[[1]] / n_to_allocate),
                        prop_n01_0 = unique(seq(prop_min[2], prop_max[2], by = prop_inc), prev_grid[[2]] / n_to_allocate),
                        prop_n10_0 = unique(seq(prop_min[3], prop_max[3], by = prop_inc), prev_grid[[3]] / n_to_allocate),
                        prop_n11_0 = unique(seq(prop_min[4], prop_max[4], by = prop_inc), prev_grid[[4]] / n_to_allocate),
                        prop_n00_1 = unique(seq(prop_min[5], prop_max[5], by = prop_inc), prev_grid[[5]] / n_to_allocate),
                        prop_n01_1 = unique(seq(prop_min[6], prop_max[6], by = prop_inc), prev_grid[[6]] / n_to_allocate),
                        prop_n10_1 = unique(seq(prop_min[7], prop_max[7], by = prop_inc), prev_grid[[7]] / n_to_allocate),
                        prop_n11_1 = unique(seq(prop_min[8], prop_max[8], by = prop_inc), prev_grid[[8]] / n_to_allocate)
                        )
    #grid <- grid[rowSums(grid) <= 1, ]
    #grid <- cbind(prop_n00_0 = (1 - rowSums(grid)), grid)

    grid <- grid[rowSums(grid) == 1, ]

    grid$n00_0 <- grid$prop_n00_0 * n_to_allocate
    grid$n01_0 <- grid$prop_n01_0 * n_to_allocate
    grid$n10_0 <- grid$prop_n10_0 * n_to_allocate
    grid$n11_0 <- grid$prop_n11_0 * n_to_allocate
    grid$n00_1 <- grid$prop_n00_1 * n_to_allocate
    grid$n01_1 <- grid$prop_n01_1 * n_to_allocate
    grid$n10_1 <- grid$prop_n10_1 * n_to_allocate
    grid$n11_1 <- grid$prop_n11_1 * n_to_allocate
  }

  grid[, ncol(grid)/2 + 1:num_strat] <- grid[, ncol(grid)/2 + 1:num_strat] + min_n

  # Make sure n = phII
  grid[, "n"] <- round(rowSums(grid[, ncol(grid)/2 + 1:num_strat]))
  grid <- grid[grid[, "n"] == phII, ]

  if (length(phI_strat) == 2) {
    # Make sure n < N in each stratum
    grid <- grid[grid[, "n0"] <= phI_strat$N0, ]
    grid <- grid[grid[, "n1"] <= phI_strat$N1, ]
    # Create columns with the P(Vi=1|sample_on)
    grid[, "pi0"] <- grid[, "n0"] / phI_strat$N0
    grid[, "pi1"] <- grid[, "n1"] / phI_strat$N1
  } else if (length(phI_strat) == 4) {
    # Make sure n < N in each stratum
    grid <- grid[grid[, "n00"] <= phI_strat$N00, ]
    grid <- grid[grid[, "n01"] <= phI_strat$N01, ]
    grid <- grid[grid[, "n10"] <= phI_strat$N10, ]
    grid <- grid[grid[, "n11"] <= phI_strat$N11, ]
    # Create columns with the P(Vi=1|sample_on)
    grid[, "pi00"] <- grid[, "n00"] / phI_strat$N00
    grid[, "pi01"] <- grid[, "n01"] / phI_strat$N01
    grid[, "pi10"] <- grid[, "n10"] / phI_strat$N10
    grid[, "pi11"] <- grid[, "n11"] / phI_strat$N11
  } else if (length(phI_strat) == 8) {
    # Make sure n < N in each stratum
    grid <- grid[grid[, "n00_0"] <= phI_strat$N00_0, ]
    grid <- grid[grid[, "n01_0"] <= phI_strat$N01_0, ]
    grid <- grid[grid[, "n10_0"] <= phI_strat$N10_0, ]
    grid <- grid[grid[, "n11_0"] <= phI_strat$N11_0, ]
    grid <- grid[grid[, "n00_1"] <= phI_strat$N00_1, ]
    grid <- grid[grid[, "n01_1"] <= phI_strat$N01_1, ]
    grid <- grid[grid[, "n10_1"] <= phI_strat$N10_1, ]
    grid <- grid[grid[, "n11_1"] <= phI_strat$N11_1, ]
    # Create columns with the P(Vi=1|sample_on)
    grid[, "pi00_0"] <- grid[, "n00_0"] / phI_strat$N00_0
    grid[, "pi01_0"] <- grid[, "n01_0"] / phI_strat$N01_0
    grid[, "pi10_0"] <- grid[, "n10_0"] / phI_strat$N10_0
    grid[, "pi11_0"] <- grid[, "n11_0"] / phI_strat$N11_0
    grid[, "pi00_1"] <- grid[, "n00_1"] / phI_strat$N00_1
    grid[, "pi01_1"] <- grid[, "n01_1"] / phI_strat$N01_1
    grid[, "pi10_1"] <- grid[, "n10_1"] / phI_strat$N10_1
    grid[, "pi11_1"] <- grid[, "n11_1"] / phI_strat$N11_1
  }
  # Remove "prop_" columns before returning grid
  grid <- grid[, -grep(pattern = "prop", x = colnames(grid))]
  grid[, "n"] <- NULL
  return(grid)
}
