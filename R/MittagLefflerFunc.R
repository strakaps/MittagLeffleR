# Copyright (c) 2015, Roberto Garrappa 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are 
# met:
#   
# * Redistributions of source code must retain the above copyright 
#   notice, this list of conditions and the following disclaimer. 
# * Redistributions in binary form must reproduce the above copyright 
#   notice, this list of conditions and the following disclaimer in 
#   the documentation and/or other materials provided with the distribution 
# * Neither the name of the Department of Mathematics - University of Bari - 
#   Italy nor the names of its contributors may be used to endorse or promote
#   products derived from this software without specific prior written 
#   permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.

# The code below is adapted from: 
# https://au.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function

#' Mittag-Leffler Function.
#' 
#' The generalized (two-parameter) Mittag-Leffer function is defined by the
#' power series
#'     \deqn{E_{\alpha,\beta} (z) = \sum_{k=0}^\infty  z^k / \Gamma(\alpha
#'     k + \beta) }
#' for complex \eqn{z} and complex \eqn{\alpha, \beta} with
#' \eqn{Real(\alpha) > 0} (only implemented for real valued parameters).
#' 
#' @rdname mlf
#' @export
#' @param z The argument (real-valued)
#' @param a,b,g Parameters of the Mittag-Leffler distribution; see Garrappa
#' @return \code{mlf} returns the value of the Mittag-Leffler function.
#' @examples
#' mlf(2,0.7)
#' @references
#' Garrappa, R. (2015). Numerical Evaluation of Two and Three Parameter
#' Mittag-Leffler Functions. SIAM Journal on Numerical Analysis, 53(3),
#'  1350–1369. \url{http://doi.org/10.1137/140971191}
#'
#' The Mittag-Leffler function. MathWorks File Exchange.
#' \url{https://au.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function}
#' 


mlf <- function(z,a,b=1,g=1) {
  #check number of inputs and handle errors
  if (nargs() < 2) {
    #handle error
  }
  
  #check inputs values and handle errors
  
  #target precision
  log_epsilon <- log(10^(-15))
  
  #LT inversion for each element of z
  E <- numeric(length(z))
  for (i in 1:length(z)) {
    if (abs(z[i]) < 10^(-15) ) {
      E[i] <- 1/gamma(b)
    }
    else {
      E[i] <- LTinversion(1,z[i],a,b,g,log_epsilon)
    }
  }
  
  return(E)
}

LTinversion <- function(t, lambda, a, b, g, log_epsilon) {
  theta <- Arg(lambda)
  kmin <- ceiling(-a/2 - theta/2/pi)
  kmax <- floor(a/2 - theta/2/pi)
  if (kmin > kmax) {
    k_vett <- c()
  }
  else {
    k_vett <- kmin:kmax
  }
  s_star <- abs(lambda)^(1/a) * exp(1i*(theta+2*k_vett*pi)/a)
  
  phi_s_star <- (Re(s_star)+abs(s_star))/2
  
  phi_s_star <- sort(phi_s_star, index.return=TRUE)$x
  index_s_star <- sort(phi_s_star, index.return=TRUE)$ix
  s_star <- s_star[index_s_star]
  
  index_save <- phi_s_star > 10^(-15)
  s_star <- s_star[index_save]
  phi_s_star <- phi_s_star[index_save]
  
  s_star <- c(0, s_star)
  phi_s_star <- c(0, phi_s_star)
  J1 <- length(s_star)
  J <- J1 - 1
  
  p <- c( max(0,-2*(a*g-b+1)) , rep(1,J)*g )
  q <- c( rep(1,J)*g , +Inf)
  phi_s_star <- c(phi_s_star, +Inf)
  
  admissible_regions <- which( (phi_s_star[1:length(phi_s_star)-1] < (log_epsilon - log(.Machine$double.eps))/t) & (phi_s_star[1:length(phi_s_star)-1] < phi_s_star[2:length(phi_s_star)]))
  
  JJ1 <- admissible_regions[length(admissible_regions)]
  mu_vett <- rep(1,JJ1)*Inf
  N_vett <- rep(1,JJ1)*Inf
  h_vett <- rep(1,JJ1)*Inf
  
  find_region <- 0
  while (!find_region) {
    for (j1 in admissible_regions) {
      if (j1 < J1) {
        temp2 <- OptimalParam_RB(t,phi_s_star[j1],phi_s_star[j1+1],p[j1],q[j1],log_epsilon)
        muj <- temp2[1]
        hj <- temp2[2]
        Nj <- temp2[3]
      }
      else {
        temp3 <- OptimalParam_RU(t,phi_s_star[j1],p[j1],log_epsilon)
        muj <- temp3[1]
        hj <- temp3[2]
        Nj <- temp3[3]
      }
      mu_vett[j1] <- muj
      h_vett[j1] <- hj
      N_vett[j1] <- Nj
    }
    if (min(N_vett) > 200) {
      log_epsilon <- log_epsilon +log(10)
    }
    else {
      find_region <- 1
    }
  }
  
  N <- min(N_vett)
  iN <- which.min(N_vett)
  mu <- mu_vett[iN]
  h <- h_vett[iN]
  
  k <- -N:N
  u <- h*k
  z <- mu*(1i*u+1)^2
  zd <- -2*mu*u + 2*mu*1i
  zexp <- exp(z*t)
  F <- z^(a*g-b)/(z^a - lambda)^g*zd
  S <- zexp*F
  Integral <- h*sum(S)/2/pi/1i
  
  if ((iN+1) > length(s_star)) {
    ss_star <- c()
  }
  else {
    ss_star <- s_star[(iN+1):length(s_star)]
  }
  Residues <- sum(1/a*(ss_star)^(1-b)*exp(t*ss_star))
  
  E <- Integral + Residues
  if (!is.complex(lambda)) {
    E <- Re(E)
  }
  return(E)
}

OptimalParam_RB <- function(t, phi_s_star_j, phi_s_star_j1, pj, qj, log_epsilon) {
  log_eps <- -36.043653389117154 #log(eps)
  fac <- 1.01
  conservative_error_analysis <- 0
  
  f_max <- exp(log_epsilon - log_eps)
  
  sq_phi_star_j <- sqrt(phi_s_star_j)
  threshold <- 2*sqrt((log_epsilon - log_eps)/t) ;
  sq_phi_star_j1 <- min(sqrt(phi_s_star_j1), threshold - sq_phi_star_j) ;
  
  if ((pj < 10^(-14)) && (qj < 10^(-14))) {
    sq_phibar_star_j <- sq_phi_star_j
    sq_phibar_star_j1 <- sq_phi_star_j1
    adm_region <- 1
  }
  
  if ((pj < 10^(-14)) && (qj >= 10^(-14))) {
    sq_phibar_star_j <- sq_phi_star_j
    if (sq_phi_star_j > 0) {
      f_min <- fac*(sq_phi_star_j/(sq_phi_star_j1-sq_phi_star_j))^qj
    }
    else {
      f_min <- fac
    }
    if (f_min < f_max) {
      f_bar <- f_min + f_min/f_max*(f_max-f_min)
      fq <- f_bar^(-1/qj)
      sq_phibar_star_j1 <- (2*sq_phi_star_j1 - fq*sq_phi_star_j)/(2+fq)
      adm_region <- 1
    }
    else {
      adm_region <- 0
    }
  }
  
  if ((pj >= 10^(-14)) && (qj < 10^(-14))) {
    sq_phibar_star_j1 <- sq_phi_star_j1
    f_min = fac*(sq_phi_star_j1/(sq_phi_star_j1 - sq_phi_star_j))^pj
    if (f_min < f_max) {
      f_bar <- f_min + f_min/f_max*(f_max - f_min)
      fp <- f_bar^(-1/pj)
      sq_phibar_star_j <- (2*sq_phi_star_j+fp*sq_phi_star_j1)/(2-fp)
      adm_region <- 1
    }
    else {
      adm_region <- 0
    }
  }
  
  if ((pj >= 10^(-14)) && (qj >= 10^(-14))) {
    f_min = fac*(sq_phi_star_j+sq_phi_star_j1)/(sq_phi_star_j1-sq_phi_star_j)^max(pj,qj)
    if (f_min < f_max) {
      f_min <- max(f_min,1.5)
      f_bar <- f_min + f_min/f_max*(f_max-f_min)
      fp <- f_bar^(-1/pj)
      fq <- f_bar^(-1/qj)
      if (!conservative_error_analysis) {
        w <- -phi_s_star_j1*t/log_epsilon
      }
      else {
        w <- -2*phi_s_star_j1*t/(log_epsilon-phi_s_star_j1*t)
      }
      den <- 2+w - (1+w)*fp + fq
      sq_phibar_star_j <- ((2+w+fq)*sq_phi_star_j + fp*sq_phi_star_j1)/den
      sq_phibar_star_j1 <- (-(1+w)*fq*sq_phi_star_j + (2+w-(1+w)*fp)*sq_phi_star_j1)/den
      adm_region <- 1
    }
    else {
      adm_region <- 0
    }
  }
  
  if (adm_region) {
    log_epsilon <- log_epsilon  - log(f_bar)
    if (!conservative_error_analysis) {
      w <- -sq_phibar_star_j1^2*t/log_epsilon
    }
    else {
      w = -2*sq_phibar_star_j1^2*t/(log_epsilon-sq_phibar_star_j1^2*t)
    }
    muj <- (((1+w)*sq_phibar_star_j + sq_phibar_star_j1)/(2+w))^2
    hj <- -2*pi/log_epsilon*(sq_phibar_star_j1-sq_phibar_star_j)/((1+w)*sq_phibar_star_j + sq_phibar_star_j1)
    Nj <- ceiling(sqrt(1-log_epsilon/t/muj)/hj)
  }
  else {
    muj <- 0
    hj <- 0
    Nj <- +Inf
  }
  
  out <- c(muj, hj, Nj)
  return(out)
  
}

OptimalParam_RU <- function(t, phi_s_star_j, pj, log_epsilon) {
  sq_phi_s_star_j <- sqrt(phi_s_star_j)
  if (phi_s_star_j > 0) {
    phibar_star_j <- phi_s_star_j*1.01
  }
  else {
    phibar_star_j <- 0.01
  }
  sq_phibar_star_j <- sqrt(phibar_star_j)
  
  f_min <- 1
  f_max <- 10
  f_tar <- 5
  
  stopp <- 0
  while (!stopp) {
    phi_t <- phibar_star_j*t
    log_eps_phi_t <- log_epsilon/phi_t
    Nj <- ceiling(phi_t/pi*(1 - 3*log_eps_phi_t/2 + sqrt(1-2*log_eps_phi_t)))
    A <- pi*Nj/phi_t
    sq_muj <- sq_phibar_star_j*abs(4-A)/abs(7-sqrt(1+12*A))
    fbar <- ((sq_phibar_star_j-sq_phi_s_star_j)/sq_muj)^(-pj)
    stopp <- (pj < 10^(-14)) || (f_min < fbar && fbar < f_max)
    if (!stopp) {
      sq_phibar_star_j <- f_tar^(-1/pj)*sq_muj + sq_phi_s_star_j
      phibar_star_j <- sq_phibar_star_j^2
    }
  }
  muj <- sq_muj^2
  hj <- (-3*A - 2 + 2*sqrt(1+12*A))/(4-A)/Nj
  
  log_eps <- log(.Machine$double.eps)
  threshold <- (log_epsilon - log_eps)/t
  if (muj > threshold) {
    if (abs(pj) < 10^(-14)) {
      Q <- 0
    }
    else {
      Q <- f_tar^(-1/pj)*sqrt(muj)
    }
    phibar_star_j <- (Q + sqrt(phi_s_star_j))^2
    if (phibar_star_j < threshold) {
      w <- sqrt(log_eps/(log_eps-log_epsilon))
      u <- sqrt(-phibar_star_j*t/log_eps)
      muj <- threshold
      Nj <- ceiling(w*log_epsilon/2/pi/(u*w-1))
      hj <- sqrt(log_eps/(log_eps - log_epsilon))/Nj
    }
    else {
      Nj <- +Inf
      hj = 0
    }
  }
  
  out <- c(muj, hj, Nj)
  return(out)
}

