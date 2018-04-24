#' Distribution functions and random number generation.
#'
#' Probability density, cumulative distribution
#' function, quantile function and random variate generation for the
#' two types of Mittag-Leffler distribution.
#' The Laplace inversion algorithm by Garrappa is used for the pdf and 
#' cdf (see 
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function}).
#'
#' The Mittag-Leffler function \code{\link{mlf}} defines two types of 
#' probability distributions:
#'
#' The \strong{first type} of Mittag-Leffler distribution assumes the Mittag-Leffler
#' function as its tail function, so that the CDF is given by
#' \deqn{F(q; \alpha, \tau) = 1 - E_{\alpha,1} (-(q/\tau)^\alpha)}
#' for \eqn{q \ge 0}, tail parameter \eqn{0 < \alpha \le 1},
#' and scale parameter \eqn{\tau > 0}.
#' Its PDF is given by
#' \deqn{f(x; \alpha, \tau) = x^{\alpha - 1} 
#' E_{\alpha,\alpha} [-(x/\tau)^\alpha] / \tau^\alpha.}
#' As \eqn{\alpha} approaches 1 from below, the Mittag-Leffler converges
#' (weakly) to the exponential
#' distribution. For \eqn{0 < \alpha < 1}, it is (very) heavy-tailed, i.e.
#' has infinite mean.
#'
#' The \strong{second type} of Mittag-Leffler distribution is defined via the
#' Laplace transform of its density f:
#' \deqn{\int_0^\infty \exp(-sx) f(x; \alpha, 1) dx = E_{\alpha,1}(-s)}
#' It is light-tailed, i.e. all its moments are finite.
#' At scale \eqn{\tau}, its density is 
#' \deqn{f(x; \alpha, \tau) = f(x/\tau; \alpha, 1) / \tau.}
#'
#' @return \code{dml} returns the density,
#'         \code{pml} returns the distribution function,
#'         \code{qml} returns the quantile function, and
#'         \code{rml} generates random variables.

#' @references
#' Haubold, H. J., Mathai, A. M., & Saxena, R. K. (2011). Mittag-Leffler
#' Functions and Their Applications. Journal of Applied Mathematics, 2011, 
#' 1–51. \url{http://doi.org/10.1155/2011/298628}
#' 
#' Mittag-Leffler distribution. (2017, May 3).
#' In Wikipedia, The Free Encyclopedia.
#' \url{https://en.wikipedia.org/w/index.php?title=Mittag-Leffler_distribution&oldid=778429885}
#' 
#' @name MLdistribution


#' @param x,q vector of quantiles.
#' @param tail tail parameter.
#' @param scale scale parameter.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param second.type logical; if FALSE (default), 
#'        first type of Mittag-Leffler distribution is assumed.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]}
#'        otherwise, \eqn{P[X > x]}
#' @rdname MittagLeffleR
#' @examples
#' dml(1, 0.8)
#' dml(1, 0.6, second.type=TRUE)
#' @name Mittag-Leffler
#' @family Mittag Leffler Distribution
#' @export
dml <- function(x,tail,scale=1,log=FALSE, second.type=FALSE){
  if (length(tail) > 1){
    stop("length(tail) must be 1.")
  }
  if (second.type==FALSE) {
    return(dml1(x,tail,scale,log))
  } else {
    y <- dml2(x/scale,tail)/scale
    if (!log) return(y) else return(log(y))
  }
}

# first type
dml1 <- function(t,tail,scale=1,log=FALSE) {
  ml <- (t^(tail-1)/(scale^tail))*mlf(-(t/scale)^tail, tail, tail, 1)
  if (log) {
    ml <- log(ml)
  }
  return(ml)
}

# second type, unit scale
dml2 <- function(u,tail) {
  # find the distribution of E(1), where E() is the inverse stable
  # subordinator; Meerschaert & Straka, Eq.(9).
  # The scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  
  h <- (1/tail)*u^(-1-1/tail)*stabledist::dstable(x = u^(-1/tail), 
        alpha = tail, beta = 1, gamma = gamma, delta = 0, pm = 1)
}


#' @examples
#' pml(2, 0.7, 1.5)
#' @name pml
#' @rdname MittagLeffleR
#' @family Mittag Leffler Distribution
#' @export
pml <- function(q, tail, scale=1, second.type=FALSE, lower.tail=TRUE, 
                log.p=FALSE) {
  # rescale
  q <- q/scale
  if (!second.type){
    p <- pml1(q,tail)
  } else {
    p <- pml2(q,tail)
  }
  if (!lower.tail) {
    p <- 1-p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

# type 1 with unit scale
pml1 <- function(q,tail) {
  p <- 1 - mlf(-(q)^tail,tail,1,1)
}

# type 2 with unit scale
pml2 <- function(q,tail) {
  # the scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  p <- stabledist::pstable(q^(-1/tail), alpha=tail, beta=1, gamma=gamma, delta=0,
                           pm=1, lower.tail = FALSE)
}

#' @rdname MittagLeffleR
#' @examples
#' qml(p = c(0.25, 0.5, 0.75), tail = 0.6, scale = 100)
#' @family Mittag Leffler Distribution
#' @export
#' @param p vector of probabilities.

qml <- function(p, tail, scale=1, second.type=FALSE, lower.tail=TRUE,
                log.p=FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1-p
  }
  if (!second.type) {
    q <- qml1(p, tail, scale)
  } else {
    q <- scale * qml2(p, tail)
  }
  return(q)
}

qml1 <- function(p,tail,scale=1) {
  x <- numeric(length(p))
  for (i in 1:length(p)) {
    qml_p <- function(t) {pml(t,tail,scale) - p[i]}
    x[i] <- stats::uniroot(qml_p, interval = c(10^-14,100),
                           extendInt="upX", tol = 1e-14)$root
  }
  return(x)
}

# type 2 with unit scale
qml2 <- function(p, tail){
  # the scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  q <- stabledist::qstable(p, alpha=tail, beta=1, gamma=gamma, delta=0, pm=1,
                           lower.tail=FALSE)^(-tail)
}

#' @rdname MittagLeffleR
#' @examples
#' rml(10, 0.7, 1)
#'

#' @export
#' @family Mittag Leffler Distribution
#' @param n number of observations. If length(n) > 1, the length is taken
#'        to be the number required.
rml <- function(n,tail,scale=1, second.type=FALSE){
  if (length(n) > 1){
    n <- length(n)
  }
  if (!second.type){
    x <- scale * rml1(n,tail)
  } else {
    x <- scale * rml2(n,tail)
  }
  return(x)
}

# unit scale; see e.g. Haubold, Mathai & Saxena (2011)
rml1 <- function(n, tail){
  # the scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  y <- stabledist::rstable(n, alpha=tail, beta=1, gamma=gamma, delta=0, pm=1)
  x <- stats::rexp(n)
  y * x^(1/tail)
}

rml2 <- function(n, tail){
  # the scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  stabledist::rstable(n, alpha=tail, beta=1, gamma=gamma, delta=0, pm=1)^(-tail)
}
