#' (Rewrite) Adaptive grid search optimization to find the optimal design (optMLE)
#' @name optMLE_grid2
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param phI_strat Phase I stratum sample sizes as a named list, dataframe, or vector
#' @param phIIa_strat For multi-wave designs, Phase II(a) stratum sample sizes as a named list, dataframe, or vector Default is \code{NULL}.
#' @param min_n Minimum stratum size to be sampled.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). All variables must be coded as binary.
#' @param indiv_score Matrix of score vectors for all parameters (returned from the \code{score} function.
#' @param return_full_grid If \code{TRUE}, all audits from all iterations of the grid search will be return. Default is \code{FALSE}.
#' @param first_step (Optional) Starting grid scale.
#' @param max_grid_size (Optional) Integer maxium for the largest starting grid allowed if \code{first_step = NULL}.
#' @return
#' \item{all_opt}{Optimal designs chosen in each iteration of the grid search.}
#' \item{min_var}{Value of the variance achieved by the optimal design in the last iteration.}
#' \item{min_var_design}{Optimal design in the last iteration.}
#' \item{findOptimal}{\code{TRUE}/\code{FALSE} for whether final iteration found an optimal design.}
#' \item{full_grid_search}{If \code{return_full_grid} = TRUE, a dataframe containing all grids from ann iterations.}
#' \item{message}{Result of the grid search, options include \code{"No valid grids"}, \code{"Singular information"}, \code{"Tie for minimum"}, \code{"Grid completed without finding minimum"}, \code{"Grid search successful"}.}
#' @export
optMLE_grid2 <- function(phI, phII, phI_strat, phIIa_strat = NULL, min_n, sample_on, indiv_score, return_full_grid = FALSE, first_step = NULL, max_grid_size = 10000) {
  # Since each of the strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - num_strat x min_n
  num_strat <- length(phI_strat)
  phi <- phII - sum(unlist(pmin(phI_strat, min_n)))

  # Create stratum IDs for merging with sampling probabilities
  indiv_score$strat <- "N"
  for (v in sample_on) {
    indiv_score$strat <- paste0(indiv_score$strat, indiv_score[, v])
  }

  # Create a dataframe to store all grid's optimal designs
  all_opt_des <- data.frame()

  it <- 1

  if (is.null(first_step)) {
    prev_step <- suggest_step(phII = phII, phI_strat = phI_strat, min_n = min_n, num_strat = num_strat, prev_grid_des = NULL, prev_delta = NULL, max_grid_size = max_grid_size) #phi
  } else {
    prev_step <- first_step
  }

  # Run initial grid search
  grid <- build_grid(delta = prev_step, phi = phi, num_strat = num_strat, phI_strat = phI_strat, phIIa_strat = phIIa_strat, min_n = min_n, prev_grid_des = NULL, prev_delta = NULL)
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

  # Save previous grid & step size
  prev_grid <- grid

  # Check for a clear minimum
  min_var_design <- prev_min <- grid[grid$Vbeta == min_var, ]
  all_opt_des <- rbind(all_opt_des,
                       cbind(grid = 1, audit_step = prev_step, min_var_design, grid_size = nrow(grid)))

  # Create dataframe to store all designs
  if (return_full_grid) { all_grids <- cbind(grid = 1, grid) }

  # Keep choosing new steps until get to 1-person scale
  keepSearching <- TRUE
  while (keepSearching) {
    it <- it + 1

    # Reset indicators
    findOptimal <- findFinalOptimal <- FALSE

    # To choose min/max for next grid, use previous grid's "optimal" design
    prev_grid_des <- prev_min[, 1:num_strat] - unlist(pmin(phI_strat, min_n))

    # If multi-wave design, and not initial grid
    ## subtract the stratum sample sizes from Phase II(a) from min_var_design
    if (!is.null(phIIa_strat)) {
      for (c in colnames(prev_grid_des)) {
        prev_grid_des[, c] <- prev_grid_des[, c] - as.numeric(phIIa_strat[c])
      }
    }

    if (max(prev_step / 2, 1) == 1) {
      prev_grid_des <- round(prev_grid_des)
    }
    # Run grid search
    grid <- build_grid(delta = max(prev_step / 2, 1), phi = phi, num_strat = num_strat, phI_strat = phI_strat, phIIa_strat = phIIa_strat, min_n = min_n, prev_grid_des = prev_grid_des, prev_delta = max(prev_step / 2, 1))
    grid$Vbeta <- apply(X = grid, MARGIN = 1, FUN = var_by_row, phI = phI, indiv_score = indiv_score, sample_on = sample_on)

    # Combine with the previous grid
    if (max(prev_step / 2, 1) > 1) {
      grid <- unique(rbind(grid, prev_grid))
    }
    min_var <- min(grid$Vbeta)

    if (is.na(min_var)) {
      all_opt_des$grid <- 1:nrow(all_opt_des)
      #all_opt_des$audit_step <- audit_steps[-length(audit_steps)]
      return(list("all_opt" = all_opt_des,
                  "min_var" = 9999,
                  "min_var_design" = NA,
                  "findOptimal" = FALSE,
                  "full_grid_search" = NA,
                  "message" = "Singular information"))
    } else if (sum(grid$Vbeta == min_var) > 1) {
      all_opt_des$grid <- 1:nrow(all_opt_des)
      #all_opt_des$audit_step <- audit_steps[-length(audit_steps)]
      return(list("all_opt" = all_opt_des,
                  "min_var" = 9999,
                  "min_var_design" = NA,
                  "findOptimal" = FALSE,
                  "full_grid_search" = NA,
                  "message" = "Tie for minimum"))
    }

    findOptimal <- sum(grid$Vbeta == min_var) == 1# & (min_var - all_opt_des$Vbeta[nrow(all_opt_des)]) < 1E-8
    # if (findOptimal &
    #     rowSums(min_var_design) > phII) { findOptimal <- TRUE }

    if (findOptimal) {
      min_var_design <- prev_min <- grid[grid$Vbeta == min_var, ]
      all_opt_des <- rbind(all_opt_des,
                           cbind(grid = it, audit_step = max(prev_step / 2, 1), min_var_design, grid_size = nrow(grid)))
      if (return_full_grid) { all_grids <- rbind(all_grids, cbind(grid = it, grid)) }
      keepSearching <- max(prev_step / 2, 1) > 1
    }

    prev_step <- prev_step / 2
    prev_grid <- grid
  }

  if (findOptimal & all_opt_des$audit_step[nrow(all_opt_des)] == 1) {findFinalOptimal <- TRUE}

  #all_opt_des$audit_step <- audit_steps

  if (!findFinalOptimal) { min_var_design[1, ] <- NA }

  if (!return_full_grid) { all_grids = NA }

  return(list("all_opt" = all_opt_des,
              "min_var" = min_var,
              "min_var_design" = min_var_design,
              "findOptimal" = findFinalOptimal,
              "full_grid_search" = all_grids,
              "message" = ifelse(findFinalOptimal, "Grid search successful", "Grid completed without finding minimum")))
}
