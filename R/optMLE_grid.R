#' Adaptive grid search optimization to find the optimal design (optMLE)
#' @name optMLE_grid
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param phI_strat Phase I stratum sample sizes as a named list, dataframe, or vector
#' @param phIIa_strat For multi-wave designs, Phase II(a) stratum sample sizes as a named list, dataframe, or vector Default is \code{NULL}.
#' @param min_n Minimum stratum size to be sampled.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 binary variables can be accommodated.
#' @param steps (Optional) Vector of step sizes for the grid search. Default is \code{NULL}, in which case the internal \code{suggest_step} function is used to choose grid scale.
#' @param indiv_score Matrix of score vectors for all parameters.
#' @param return_full_grid If \code{TRUE}, all audits from all iterations of the grid search will be return. Default is \code{FALSE}.
#' @param max_grid_size Integer maxium for the largest grids that will be searched.
#' @return
#' \item{all_opt}{Optimal designs chosen in each iteration of the grid search.}
#' \item{min_var}{Value of the variance achieved by the optimal design in the last iteration.}
#' \item{min_var_design}{Optimal design in the last iteration.}
#' \item{findOptimal}{\code{TRUE}/\code{FALSE} for whether final iteration found an optimal design.}
#' \item{full_grid_search}{If \code{return_full_grid} = TRUE, a dataframe containing all grids from ann iterations.}
#' \item{message}{Result of the grid search, options include \code{"No valid grids"}, \code{"Singular information"}, \code{"Tie for minimum"}, \code{"Grid completed without finding minimum"}, \code{"Grid search successful"}.}
#' @export
optMLE_grid <- function(phI, phII, phI_strat, phIIa_strat = NULL, min_n, sample_on, steps = NULL, indiv_score, return_full_grid = FALSE, max_grid_size = 10000) {
  # Since each of the strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - num_strat x min_n
  num_strat <- 2 ^ length(sample_on)
  phi <- phII - num_strat * min_n

  if (is.null(steps)) {
    # Initial audit step size
    audit_steps <- as.vector(suggest_step(phII = phII, phI_strat = phI_strat, min_n = min_n, num_strat = num_strat, prev_grid_des = NULL, prev_delta = NULL, max_grid_size = max_grid_size))
    if (audit_steps == Inf) {
      return(list("all_opt" = NA,
                  "min_var" = 9999,
                  "min_var_design" = NA,
                  "findOptimal" = FALSE,
                  "full_grid_search" = NA,
                  "message" = "No valid grids"))
    }
  } else {
    audit_steps <- steps
  }

  # Create a dataframe to store all grid's optimal designs
  all_opt_des <- data.frame()

  if (!(num_strat %in% c(2, 4, 8))) {
    return(warning("optMLE_grid() handles sampling on 2, 4, or 8 strata."))
  }

  it <- 1

  # Run initial grid search
  grid <- build_grid(delta = audit_steps[it], phi = phi, num_strat = num_strat, phI_strat = phI_strat, phIIa_strat = phIIa_strat, min_n = min_n, prev_grid_des = NULL, prev_delta = NULL)

  if (any(grid[, grep("pi", colnames(grid))] > 1)) {
    return(warning("Invalid grid values - sampling more than available."))
  }

  grid$Vbeta <- apply(X = grid, MARGIN = 1, FUN = var_by_row, phI = phI, indiv_score = indiv_score, sample_on = sample_on)
  min_var <- min(grid$Vbeta)

  if (is.na(min_var)) {
    return(list("all_opt" = NA,
                "min_var" = 9999,
                "min_var_design" = NA,
                "findOptimal" = FALSE,
                "full_grid_search" = NA,
                "message" = "Singular information"))
  } else if (sum(grid$Vbeta == min_var) > 1) {
    return(list("all_opt" = all_opt_des,
                "min_var" = 9999,
                "min_var_design" = NA,
                "findOptimal" = FALSE,
                "full_grid_search" = NA,
                "message" = "Tie for minimum"))
  }

  # Check for a clear minimum
  min_var_design <- prev_min <- grid[grid$Vbeta == min_var, ]
  all_opt_des <- rbind(all_opt_des,
                       cbind(grid = 1, audit_step = NA, min_var_design, grid_size = nrow(grid)))

  # Create dataframe to store all designs
  if (return_full_grid) { all_grids <- cbind(grid = 1, grid) }

  # If first grid is at 1-person scale, check for optimal
  if (audit_steps[it] == 1 & sum(grid$Vbeta == min_var) == 1) {
    findFinalOptimal <- TRUE
  }

  # Keep choosing new steps until get to 1-person scale
  keepSearching <- ifelse(is.null(steps), audit_steps[it] > 1, it < length(audit_steps))
  while (keepSearching) {
    it <- it + 1

    # Reset indicators
    findOptimal <- findFinalOptimal <- FALSE

    # To choose min/max for next grid, use previous grid's "optimal" design
    prev_grid_des <- prev_min[, 1:num_strat] - min_n

    # If multi-wave design, and not initial grid
    ## subtract the stratum sample sizes from Phase II(a) from min_var_design
    if (!is.null(phIIa_strat)) {
      for (c in colnames(prev_grid_des)) {
        prev_grid_des[, c] <- prev_grid_des[, c] - as.numeric(phIIa_strat[c])
      }
    }

    # Initial audit step size
    if (is.null(steps)) {
      audit_steps <- append(audit_steps,
                            suggest_step(phII = phII, phI_strat = phI_strat, min_n = min_n, num_strat = num_strat, prev_grid_des = prev_grid_des, prev_delta = audit_steps[(it - 1)], max_grid_size = max_grid_size))
    }

    if (any(audit_steps == Inf)) {
      all_opt_des$grid <- 1:nrow(all_opt_des)
      all_opt_des$audit_step <- audit_steps[-length(audit_steps)]
      return(list("all_opt" = all_opt_des,
                  "min_var" = 9999,
                  "min_var_design" = NA,
                  "findOptimal" = FALSE,
                  "full_grid_search" = NA,
                  "message" = "No valid grids"))
    }

    # Run initial grid search
    grid <- build_grid(delta = audit_steps[it], phi = phi, num_strat = num_strat, phI_strat = phI_strat, phIIa_strat = phIIa_strat, min_n = min_n, prev_grid_des = prev_grid_des, prev_delta = audit_steps[it - 1])

    if (any(grid[, grep("pi", colnames(grid))] > 1)) {
      return(warning("Invalid grid values - sampling more than available."))
    }

    grid$Vbeta <- apply(X = grid, MARGIN = 1, FUN = var_by_row, phI = phI, indiv_score = indiv_score, sample_on = sample_on)
    min_var <- min(grid$Vbeta)

    if (is.na(min_var)) {
      all_opt_des$grid <- 1:nrow(all_opt_des)
      all_opt_des$audit_step <- audit_steps[-length(audit_steps)]
      return(list("all_opt" = all_opt_des,
                  "min_var" = 9999,
                  "min_var_design" = NA,
                  "findOptimal" = FALSE,
                  "full_grid_search" = NA,
                  "message" = "Singular information"))
    } else if (sum(grid$Vbeta == min_var) > 1) {
      all_opt_des$grid <- 1:nrow(all_opt_des)
      all_opt_des$audit_step <- audit_steps[-length(audit_steps)]
      return(list("all_opt" = all_opt_des,
                  "min_var" = 9999,
                  "min_var_design" = NA,
                  "findOptimal" = FALSE,
                  "full_grid_search" = NA,
                  "message" = "Tie for minimum"))
    }

    findOptimal <- sum(grid$Vbeta == min_var) == 1 & (min_var - all_opt_des$Vbeta[nrow(all_opt_des)]) < 1E-8

    if (findOptimal) {
      min_var_design <- prev_min <- grid[grid$Vbeta == min_var, ]
      all_opt_des <- rbind(all_opt_des,
                           cbind(grid = 1, audit_step = NA, min_var_design, grid_size = nrow(grid)))
      if (return_full_grid) { all_grids <- rbind(all_grids, cbind(grid = NA, grid)) }
      if (audit_steps[length(audit_steps)] == 1) { findFinalOptimal <- TRUE }
    }
    keepSearching <- ifelse(is.null(steps), audit_steps[it] > 1, it < length(audit_steps))
  }

  all_opt_des$grid <- 1:nrow(all_opt_des)
  all_opt_des$audit_step <- audit_steps

  if (!findFinalOptimal) { min_var_design[1, ] <- NA }

  if (!return_full_grid) { all_grids = NA }

  return(list("all_opt" = all_opt_des,
              "min_var" = min_var,
              "min_var_design" = min_var_design,
              "findOptimal" = findFinalOptimal,
              "full_grid_search" = all_grids,
              "message" = ifelse(findFinalOptimal, "Grid search successful", "Grid completed without finding minimum")))
}
