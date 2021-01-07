#' Logistic regression model probability P(Y|mu)
#' @param y number or numeric vector input for value of the logistic outcome
#' @param mu value of the linear predictor
#' @name prob_logistic
#' @export
prob_logistic <- function(y, mu) {
  return(y * sigmoid(mu) + (1 - y) * (1 - sigmoid(mu)))
}
