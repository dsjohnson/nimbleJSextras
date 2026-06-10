make_GRBF <- function(x, num_knots){
  sig <- diff(range(x))/(num_knots-1)
  rbf <- function(X,Y){dnorm(X,Y,sd=sig)}
  y <- seq(min(x)-sig,max(x)+sig, length=num_knots)
  B <- outer(x,y,rbf)
  m <- colMeans(B)
  B <- sweep(B, 2, m)
  return(list(basis=B, knots=y))
}

make_honu_Gamma <- nimbleFunction(
  run = function(xi_tilde=double(), phi=double(), alpha=double(), p=double(), tau=double(1)) {
    returnType(double(2))
    m <- length(tau)
    gdim <- 6+2*m
    # submatrices

    psi_null_R <- c(1-alpha*p, alpha*p, 0, 0)

    psi_F_F <- nimMatrix(0,2*m,2*m)
    for(j in 1:(m-1)){
      psi_F_F[j,j+1] <- 1-tau[j]
      psi_F_F[m+j,m+j+1] <- 1-tau[j]
    }
    psi_F_F[m,m] <- 1-tau[m]
    psi_F_F[2*m,2*m] <- 1-tau[m]

    psi_F_R <- nimMatrix(0,2*m, 4)
    for(j in 1:m){
      psi_F_R[j,1] <- (1-alpha*p)*tau[j]
      psi_F_R[j,2] <- alpha*p*tau[j]
      psi_F_R[m+j,3] <- (1-p)*tau[j]
      psi_F_R[m+j,4] <- p*tau[j]
    }

    psi_R_F <- nimMatrix(0,4,2*m)
    psi_R_F[1,1] <- 1
    psi_R_F[2:4,m+1] <- 1

    Gamma <- nimMatrix(0,gdim, gdim)
    Gamma[1,1] <- 1-xi_tilde
    Gamma[1,2:5] <- xi_tilde*psi_null_R
    Gamma[2:5, 6:(5+2*m)] <- phi*psi_R_F
    Gamma[2:5,gdim] <- 1-phi
    Gamma[6:(gdim-1),2:5] <- phi*psi_F_R
    Gamma[6:(gdim-1), 6:(gdim-1)] <- phi*psi_F_F
    Gamma[6:(gdim-1),gdim] <- 1-phi
    Gamma[gdim,gdim] <- 1

    return(Gamma)
  })

make_honu_pi <- nimbleFunction(
  run = function(xi_tilde=double(), alpha=double(), m=double()) {
    returnType(double(1))
    pi <- rep(0,6+2*m)
    pi[1] <- 1-xi_tilde
    pi[2:5] <- xi_tilde*c(1-alpha, alpha, 0, 0)
    return(pi)
  })


make_honu_P <- nimbleFunction(
  run = function(m=integer()) {
    returnType(double(2))
    Pmat <- nimMatrix(0,nrow=6+2*m, ncol=3)
    Pmat[,1] <- 1
    Pmat[3,1] <- 0
    Pmat[5,1] <- 0
    Pmat[3,2] <- 1
    Pmat[5,3] <- 1
    return(Pmat)
  })


