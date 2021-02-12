#' Generate the sequences from window lowerbound to upperbound to fill grid
#' @name grid_vals
#' @export
grid_vals <- Vectorize(FUN = function(x, y) seq(x, y), SIMPLIFY = F)
