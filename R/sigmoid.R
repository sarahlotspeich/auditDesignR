#' Sigmoid function: f(x) = (1 + exp(x)) ^ (-1)
#' @param x Input number or numeric vector
#' @name sigmoid
#' @export
sigmoid <- function(x) {
  return((1 + exp(- x)) ^ (- 1))
}
