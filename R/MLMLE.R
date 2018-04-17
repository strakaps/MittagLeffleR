#' Maximum Likelihood Estimation of the Mittag-Leffler distribution
#'
#' Optimizes the bivariate likelihood function via 
#' \code{\link{optim}}. Uses \code{\link{logMomentEstimator}}
#' for initial parameter values.
#'
#' @param data Vector of class "numeric"
#' @param ... Additional parameters passed on to \code{\link{optim}}. 
#' @return A named vector of class "numeric"
#' @export
mlml <- function(data) {
  # theta_orig: the parameter (tail, scale)
  # theta: a transformed parameter so that optimisation is unconstrained
  theta_orig <- function(theta)
    c(1 / (1 + exp(-theta[1])), exp(theta[2]))
  theta <-
    function(theta_orig)
      c(-log(1 / theta_orig[1] - 1), log(theta_orig[2]))
  # find initial estimate
  theta_orig_init <- logMomentEstimator(data)[1:2]
  theta_init <- theta(theta_orig_init)

  # run optimization in unconstrained parameters
  log_l <- function(theta) {
    theta <- theta_orig(theta)
    - sum(dml(data, theta[1], theta[2], log = TRUE))
  }
  ml_theta <- stats::optim(theta_init, fn = log_l)$par
  tail <- 1 / (1 + exp(-ml_theta[1]))
  scale <- exp(ml_theta[2])
  c(tail, scale)
}
