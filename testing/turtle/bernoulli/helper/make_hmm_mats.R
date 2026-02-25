make_GRBF <- function(x, num_knots){
  sig <- diff(range(x))/(num_knots-1)
  rbf <- function(X,Y){dnorm(X,Y,sd=sig)}
  y <- seq(min(x),max(x), length=num_knots)
  B <- outer(x,y,rbf)
}


make_dt_pois <- nimbleFunction(
  run = function(n=integer(), rate=double(), shift=integer(default=1)){
    returnType(double(1))
    if(shift<1) stop("Poisson 'shift' must be >= 1")
    x <- 1:n
    out <- dpois(x-shift, rate)
    return(out)
  })


make_honu_Gamma <- nimbleFunction(
  run = function(xi=double(), phi=double(), p=double(1), alpha=double(), zeta=double(1)) {
    returnType(double(2))
    nd <- length(zeta)
    Fzeta <- zeta
    for(j in 2:nd){Fzeta[j] <- sum(zeta[1:j])}
    zeta_tilde <- zeta
    zeta_tilde[2:nd] <- zeta[2:nd]/(1-Fzeta[1:(nd-1)])
    gdim <- 6+2*nd
    # submatrices

    alpha_null_R <- c((1-alpha)*(1-p[1]), (1-alpha)*p[1], alpha*(1-p[2]), alpha*p[2])

    psi_F_F <- nimMatrix(0,2*nd,2*nd)
    for(j in 1:(nd-1)){
      psi_F_F[j,j+1] <- 1-zeta_tilde[j]
      psi_F_F[nd+j,nd+j+1] <- 1-zeta_tilde[j]
    }
    psi_F_F[nd,nd] <- 1-zeta_tilde[nd]
    psi_F_F[2*nd,2*nd] <- 1-zeta_tilde[nd]

    psi_F_R <- nimMatrix(0,2*nd, 4)
    for(j in 1:nd){
      psi_F_R[j,1] <- (1-p[1])*zeta_tilde[j]
      psi_F_R[j,2] <- p[1]*zeta_tilde[j]
      psi_F_R[nd+j,3] <- (1-p[2])*zeta_tilde[j]
      psi_F_R[nd+j,4] <- p[2]*zeta_tilde[j]
    }

    psi_R_F <- nimMatrix(0,4,2*nd)
    psi_R_F[1,1] <- 1
    psi_R_F[2:4,nd+1] <- 1

    Gamma <- nimMatrix(0,gdim, gdim)
    Gamma[1,1] <- 1-xi
    Gamma[1,2:5] <- xi*alpha_null_R
    Gamma[2:5, 6:(5+2*nd)] <- phi*psi_R_F
    Gamma[2:5,gdim] <- 1-phi
    Gamma[6:(gdim-1),2:5] <- phi*psi_F_R
    Gamma[6:(gdim-1), 6:(gdim-1)] <- phi*psi_F_F
    Gamma[6:(gdim-1),gdim] <- 1-phi
    Gamma[gdim,gdim] <- 1

    return(Gamma)
  })


make_honu_pi <- nimbleFunction(
  run = function(xi=double(), p=double(1), alpha=double(), nd=double()) {
    returnType(double(1))
    pi <- rep(0,6+2*nd)
    pi[1] <- 1-xi
    pi[2:5] <- xi*c((1-alpha)*(1-p[1]), (1-alpha)*p[1], alpha*(1-p[2]), alpha*p[2])
    return(pi)
  })


make_honu_P <- nimbleFunction(
  run = function(nd=integer()) {
    returnType(double(2))
    Pmat <- nimMatrix(0,nrow=6+2*nd, ncol=3)
    Pmat[,1] <- 1
    Pmat[3,1] <- 0
    Pmat[5,1] <- 0
    Pmat[3,2] <- 1
    Pmat[5,3] <- 1
    return(Pmat)
  })


