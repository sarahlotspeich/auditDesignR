library(testthat)
library(auditDesignR)

test_that("build_grid from optMLE", {
	phI = N
	phII = n
	phI_strat = stratN
	min_n = 10
	sample_on = c("Ystar", "Xstar", "Z")
	indiv_score = s
	max_grid_size = 10000
	phIIa_strat = NULL

	set.seed(918)

	num_strat <- length(phI_strat)
	phi <- phII - sum(unlist(pmin(phI_strat, min_n)))

	  # Create stratum IDs for merging with sampling probabilities
	  indiv_score$strat <- "N"
	  for (v in sample_on) {
	  	indiv_score$strat <- paste0(indiv_score$strat, indiv_score[, v])
	  }

	  # if (is.null(steps)) {
	    # Initial audit step size
	    audit_steps <- as.vector(suggest_step(phII = phII, phI_strat = phI_strat, min_n = min_n, num_strat = num_strat, prev_grid_des = NULL, prev_delta = NULL, max_grid_size = max_grid_size))
	    expect_not_equal(audit_steps, 9999)
	    # if (audit_steps == 9999) {
	    # 	return(list("all_opt" = NA,
	    # 		"min_var" = 9999,
	    # 		"min_var_design" = NA,
	    # 		"findOptimal" = FALSE,
	    # 		"full_grid_search" = NA,
	    # 		"message" = "No valid grids"))
	    # }
	    # } else {
	    # 	audit_steps <- steps
	    # }

	  # Create a dataframe to store all grid's optimal designs
	  all_opt_des <- data.frame()

	  it <- 1

	  # Run initial grid search
	  grid <- build_grid(delta = audit_steps[it], phi = phi, num_strat = num_strat, phI_strat = phI_strat, phIIa_strat = phIIa_strat, min_n = min_n, prev_grid_des = NULL, prev_delta = NULL)

	  expect_equal(dim(grid), c(6435, 16))
	  expect_equal(as.numeric(grid[1,]), c(330,10,10,10,10,10,10,10,0.109054857898215,0.00834028356964137,0.00509424350483953,0.00785545954438335,0.0137741046831956,0.0129198966408269,0.0268817204301075,0.0149925037481259))
	  })