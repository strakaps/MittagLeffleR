# From Dexter Cahoy, http://www2.latech.edu/~dcahoy/
#######################
#Random number generation and  estimation for the Mittag-Leffler distribution  
#and fractional Poisson process
#SOURCE: (1) Parameter estimation for fractional Poisson processes (with V Uchaikin and W Woyczynski). 
#Journal of Statistical Planning and Inference, 140(11), 3106-3120, Nov 2010.
#(2)Estimation of Mittag-Leffler parameters. Communications in Statistics-Simulation and Computation, 42(2), 303-315, 2013
#Updated last 02/03/2016.
#Email me at dcahoy@latech.edu if you have questions, suggestions, etc.
########################


#' Log-Moments Estimator for the Mittag-Leffler Distribution.
#' 
#' Tail and scale parameter of the Mittag-Leffler distribution are estimated
#' by matching with the first two empirical log-moments.
#' 
#' @param x     A vector of non-negative data.
#' @param alpha Confidence level is 1 - alpha.
#' @rdname logMomentEstimator
#' @export
#' @examples
#' logMomentEstimator(rml(n = 100000, scale = 0.03, tail = 0.99))
#' @return 
#' A named vector with entries (nu     delta      nuLo      nuHi   deltaLo   
#' deltaHi)

logMomentEstimator = function (x, a=0.05) {
  EULER.C = 0.57721566490153286
  log.x = log(x)
  m = mean(log.x)
  s.2 = var(log.x)
  nu = pi/sqrt(3*(s.2 + pi^2/6))
  delta = exp(m + EULER.C)
  n=length(x)
  
  se.nu=sqrt( (nu^2)*(32-20*nu^2-nu^4)/(40*n)  )
  zcv=qnorm(1-a/2,0,1)  
  l.nu= nu -zcv*se.nu
  u.nu =  nu + zcv*se.nu
  
  se.delta = sqrt(((pi^2*delta^2)/(6*n))*((2/nu^2) - 1)) #delta
  l.delta= delta -zcv*se.delta
  u.delta =  delta + zcv*se.delta
  
  return(c(nu = nu, delta = delta, nuLo = l.nu, nuHi = u.nu, 
           deltaLo = l.delta, deltaHi = u.delta))      
}

#Generate Mittag-Leffler distributed data 
#or inter-event or  inter-jump times of a fractional Poisson process
# dat=mittag.leffler(n=500,nu=0.2, mu=2)

#Point and interval estimates of  nu and mu. 
# ml.par.est(dat, a=0.05)

