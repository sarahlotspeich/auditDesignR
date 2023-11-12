#' Adaptive grid search optimization to find the optimal design (optMLE)
#' @name optMLE_grid
#' @param phI Phase I sample size.
#' @param phII Phase II sample size.
#' @param phI_strat Phase I stratum sample sizes as a named list, dataframe, or vector
#' @param phIIa_strat For multi-wave designs, Phase II(a) stratum sample sizes as a named list, dataframe, or vector. Default is \code{NULL}, assuming a single-wave design.
#' @param min_n Minimum stratum size to be sampled.
#' @param sample_on Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). All variables must be coded as binary.
#' @param steps (Optional) Vector of step sizes for the grid search. Default is \code{NULL}, in which case the internal \code{suggest_step} function is used to choose grid scale.
#' @param indiv_score Matrix of score vectors for all parameters (returned from the \code{score} function.
#' @param return_full_grid (Optional) If \code{TRUE}, all audits from all iterations of the grid search will be return. Default is \code{FALSE}.
#' @param max_grid_size Integer maximum for the largest grids that will be searched.
#' @param stop_when (Optional) Criterion used to decide when to stop the grid search. Default is \code{"step_size"}, stopping when the grid is at a 1-person level, but other option is \code{"percent_change"}, stopping when the variances of successive optimal designs change by less than 1\%.
#' @return
#' \item{all_opt}{Optimal designs chosen in each iteration of the grid search.}
#' \item{min_var}{Value of the variance achieved by the optimal design in the last iteration.}
#' \item{min_var_design}{Optimal design in the last iteration.}
#' \item{findOptimal}{\code{TRUE}/\code{FALSE} for whether final iteration found an optimal design.}
#' \item{full_grid_search}{If \code{return_full_grid} = TRUE, a dataframe containing all grids from ann iterations.}
#' \item{message}{Result of the grid search, options include \code{"No valid grids"}, \code{"Singular information"}, \code{"Tie for minimum"}, \code{"Grid completed without finding minimum"}, \code{"Grid search successful"}.}
#' @export
optMLE_grid <- function(phI, phII, phI_strat, phIIa_strat = NULL, min_n, sample_on, steps = NULL, indiv_score, return_full_grid = FALSE, max_grid_size = 10000, stop_when = "step_size") {
  # Since each of the strata must have >= min_n subjects
  ## The number that can be optimally allocated between them is only phII - num_strat x min_n
  num_strat <- length(phI_strat)
  phi <- phII - sum(unlist(pmin(phI_strat, min_n)))

  # Create stratum IDs for merging with sampling probabilities
  indiv_score$strat <- "N"
  for (v in sample_on) {
    indiv_score$strat <- paste0(indiv_score$strat, indiv_score[, v])
  }

  if (is.null(steps)) {
    # Initial audit step size
    audit_steps <- as.vector(suggest_step(phII = phII, phI_strat = phI_strat, min_n = min_n, num_strat = num_strat, prev_grid_des = NULL, prev_delta = NULL, max_grid_size = max_grid_size))
    if (audit_steps == 9999) {
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

  it <- 1

  # Run initial grid search
  grid <- build_grid(delta = audit_steps[it], phi = phi, num_strat = num_strat, phI_strat = phI_strat, phIIa_strat = phIIa_strat, min_n = min_n, prev_grid_des = NULL, prev_delta = NULL)

  if (any(grid[, grep("pi", colnames(grid))] > 1)) {
    return(warning("Invalid grid values - sampling more than available."))
  }

  grid$Vbeta <- apply(X = grid, 
                      MARGIN = 1, 
                      FUN = var_by_row, 
                      phI = phI, 
                      indiv_score = indiv_score, 
                      sample_on = sample_on)
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
  keepSearching <- ifelse(test = is.null(steps), 
                          yes = audit_steps[it] > 1, 
                          no = it < length(audit_steps))
  while (keepSearching) {
    it <- it + 1

    # Reset indicators
    findOptimal <- findFinalOptimal <- FALSE

    # To choose min/max for next grid, use previous grid's "optimal" design
    prev_grid_des <- prev_min[, 1:num_strat] - unlist(pmin(phI_strat, min_n))

    # If multi-wave design, and not initial grid,
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

    if (any(audit_steps == 9999)) {
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
    grid <- build_grid(delta = audit_steps[it], 
                       phi = phi, 
                       num_strat = num_strat, 
                       phI_strat = phI_strat, 
                       phIIa_strat = phIIa_strat, 
                       min_n = min_n, 
                       prev_grid_des = prev_grid_des, 
                       prev_delta = audit_steps[it - 1])

    if (any(grid[, grep("pi", colnames(grid))] > 1)) {
      return(warning("Invalid grid values - sampling more than available."))
    }

    grid$Vbeta <- apply(X = grid, 
                        MARGIN = 1, 
                        FUN = var_by_row, 
                        phI = phI, 
                        indiv_score = indiv_score, 
                        sample_on = sample_on)
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

    # Check for a unique "best" design from current iteration
    findOptimal <- sum(grid$Vbeta == min_var) == 1 & (min_var - all_opt_des$Vbeta[nrow(all_opt_des)]) < 1E-8
    uniqueMinimum <- sum(grid$Vbeta == min_var) == 1
    monotoneDecreasing <- (min_var - all_opt_des$Vbeta[nrow(all_opt_des)]) > 1E-8
    auditSize <- rowSums(min_var_design) > phII
    findOptimal <- uniqueMinimum & monotoneDecreasing & auditSize

    # Depending on stopping rule, decide whether to keep going or terminate search.
    if (findOptimal) {
      min_var_design <- prev_min <- grid[grid$Vbeta == min_var, ]
      all_opt_des <- rbind(all_opt_des,
                           cbind(grid = it, 
                                 audit_step = audit_steps[it], 
                                 min_var_design, 
                                 grid_size = nrow(grid)))
      if (return_full_grid) { 
        all_grids <- rbind(all_grids, 
                           cbind(grid = it, 
                                 grid)) 
      }
      
      if (stop_when == "step_size") {
        findFinalOptimal <- audit_steps[length(audit_steps)] == 1
        keepSearching <- ifelse(test = is.null(steps), 
                                yes = !findFinalOptimal & audit_steps[it] > 1, 
                                no = !findFinalOptimal & it < length(steps))  
      } else if (stop_when == "percent_change") {
        findFinalOptimal <- abs((all_opt_des$Vbeta[nrow(all_opt_des)] - min_var)) / all_opt_des$Vbeta[nrow(all_opt_des)] < 0.01
        keepSearching <- !findFinalOptimal & audit_steps[it] > 1 
      }
    }
  }

  # if (findFinalOptimal) {
  #   all_opt_des$grid <- 1:nrow(all_opt_des)
  # }
  # all_opt_des$audit_step <- audit_steps

  if (!findFinalOptimal) { min_var_design[1, ] <- NA }

  if (!return_full_grid) { all_grids = NA }

  return(list("all_opt" = all_opt_des,
              "min_var" = min_var,
              "min_var_design" = min_var_design,
              "findOptimal" = findFinalOptimal,
              "full_grid_search" = all_grids,
              "message" = ifelse(test = findFinalOptimal, 
                                 yes = "Grid search successful", 
                                 no = "Grid completed without finding minimum")))
}
