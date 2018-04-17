#' Maximum Likelihood Estimation of the Mittag-Leffler distribution
#' 
#'

mlml <- function(X) {
  log_l <- function(theta) {
    #transform parameters so can do optimization unconstrained
    theta[1] <- 1/(1+exp(-theta[1]))
    theta[2] <- exp(theta[2])
    -sum(log(dml(X,theta[1],theta[2])))
  }
  ml_theta <- stats::optim(c(0.5,0.5), fn=log_l)$par
  ml_a <- 1/(1 + exp(-ml_theta[1]))
  ml_d <- exp(ml_theta[2])
  print(paste("tail =", ml_a, "scale =", ml_d))
}
mlml(rml(n = 100, tail = 0.9, scale = 2))
