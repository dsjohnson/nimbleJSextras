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


# make_dt_Gamma <- nimbleFunction(
#   run = function(xi=double(), phi = double(0), zeta=double(1)) {
#     returnType(double(2))
#
#     nd <- length(zeta)
#     Fzeta <- zeta
#     for(j in 2:nd){Fzeta[j] <- sum(zeta[1:j])}
#     zeta_tilde <- zeta
#     zeta_tilde[2:nd] <- zeta[2:nd]/(1-Fzeta[1:(nd-1)])
#
#     Gamma <- nimMatrix(0,3+nd,3+nd)
#     Gamma[1,1] <- 1-xi
#     Gamma[1,2] <- xi
#     Gamma[2,3] <- phi
#     Gamma[3:(nd+2),2] <- phi*zeta_tilde
#     for(j in 1:(nd-1)){
#       Gamma[(2+j),(3+j)] <- phi*(1-zeta_tilde[j])
#     }
#     Gamma[2+nd,2+nd] <- phi*(1-zeta_tilde[nd])
#     Gamma[2:(2+nd), 3+nd] <- 1-phi
#     Gamma[3+nd, 3+nd] <- 1
#     return(Gamma)
#   })

make_honu_Gamma <- nimbleFunction(
  run = function(xi=double(), phi = double(0), zeta=double(1)) {
    returnType(double(2))

    nd <- length(zeta)
    Fzeta <- zeta
    for(j in 2:nd){Fzeta[j] <- sum(zeta[1:j])}
    zeta_tilde <- zeta
    zeta_tilde[2:nd] <- zeta[2:nd]/(1-Fzeta[1:(nd-1)])

    Gamma <- nimMatrix(0,4+nd,4+nd)
    Gamma[1,1] <- 1-xi
    Gamma[1,2] <- xi
    Gamma[2:3,4] <- phi
    Gamma[2:3,4+nd] <- 1-phi
    Gamma[4:(4+nd-1),3] <- phi*zeta_tilde
    for(j in 0:(nd-2)){
      Gamma[(4+j),(4+j+1)] <- phi*(1-zeta_tilde[j+1])
    }
    Gamma[3+nd,3+nd] <- phi*(1-zeta_tilde[nd])
    Gamma[4:(3+nd),4+nd] <- 1-phi
    Gamma[4+nd, 4+nd] <- 1
    return(Gamma)
  })


# make_dt_pi <- nimbleFunction(
#   run = function(xi=double(), nd=integer()) {
#     returnType(double(1))
#     pi <- rep(0,3+nd)
#     pi[1] <- 1-xi
#     pi[2] <- xi
#     return(pi)
#   })
make_honu_pi <- nimbleFunction(
  run = function(xi=double(), nd=integer()) {
    returnType(double(1))
    pi <- rep(0,4+nd)
    pi[1] <- 1-xi
    pi[2] <- xi
    return(pi)
  })

# make_dt_Pmats <- nimbleFunction(
#   run = function(p=double(), nd=integer()) {
#     returnType(double(2))
#     Pmat <- nimMatrix(nrow=nd+3, ncol=2)
#     Pmat[,1] <- 1
#     Pmat[2,1] <- 1-p
#     Pmat[2,2] <- p
#     return(Pmat)
#   })
make_honu_Rmat <- nimbleFunction(
  run = function(r=double(2), nd=integer()) {
    returnType(double(2))
    n_states <- 4+nd
    K <- dim(r)[1]
    Rmat <- nimMatrix(0, K, n_states)
    Rmat[,1] <- 0
    Rmat[,2] <- r[,1]
    Rmat[,3] <- r[,2]
    return(Rmat)
  })


