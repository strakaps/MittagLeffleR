#' Maximum Likelihood Estimation of the Mittag-Leffler distribution
#'
#' Optimizes the bivariate loglikelihood of the Mittag-Leffler distribution
#'  via \code{\link{optim}}. Uses \code{\link{logMomentEstimator}}
#' for initial parameter values.
#'
#' @param data Vector of class "numeric"
#' @param ... Additional parameters passed on to \code{\link{optim}}.
#' @return The output of \code{\link{optim}}.
#' @examples
#' library(magrittr)
#' rml(n = 100, tail = 0.8, scale = 1000) %>% mlmle()
#' @export
mlmle <- function(data, ...) {
  # theta_orig: the parameter (tail, scale)
  # theta: a transformed parameter so that optimisation is unconstrained
  theta_orig <- function(theta)
    c(1.1 / (1 + exp(-theta[1])), exp(theta[2]))
  theta <-
    function(theta_orig)
      c(-log(1.1 / theta_orig[1] - 1), log(theta_orig[2]))
  # find initial estimate
  theta_orig_init <- logMomentEstimator(data)[1:2]
  theta_orig_init[1] <- min(theta_orig_init[1], 0.99)
  theta_init <- theta(theta_orig_init)
  
  # run optimization in unconstrained parameters
  log_l <- function(theta) {
    theta <- theta_orig(theta)
    - sum(dml(data, theta[1], theta[2], log = TRUE))
  }
  opt_out <- stats::optim(theta_init, fn = log_l, ...)
  alt_loglik  <- -opt_out$value
  list(
    par = theta_orig(opt_out$par),
    loglik = alt_loglik
  )
}
