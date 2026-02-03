

load_dt_pois <- nimbleFunction(
  run = function(n=integer(), rate=double(), shift=integer(default=1)){
    returnType(double(1))
    if(shift<1) stop("Poisson 'shift' must be >= 1")
    x <- 1:n
    out <- dpois(x-shift, rate)
    return(out)
  })


load_dt_Gamma <- nimbleFunction(
  run = function(xi=double(), phi = double(0), zeta=double(1)) {
    returnType(double(2))

    nd <- length(zeta)
    Fzeta <- zeta
    for(j in 2:nd){Fzeta[j] <- sum(zeta[1:j])}
    zeta_tilde <- zeta
    zeta_tilde[2:nd] <- zeta[2:nd]/(1-Fzeta[1:(nd-1)])

    Gamma <- nimMatrix(0,3+nd,3+nd)
    Gamma[1,1] <- 1-xi
    Gamma[1,2] <- xi
    Gamma[2,3] <- phi
    Gamma[3:(nd+2),2] <- phi*zeta_tilde
    for(j in 1:(nd-1)){
      Gamma[(2+j),(3+j)] <- phi*(1-zeta_tilde[j])
    }
    Gamma[2+nd,2+nd] <- phi*(1-zeta_tilde[nd])
    Gamma[2:(2+nd), 3+nd] <- 1-phi
    Gamma[3+nd, 3+nd] <- 1
    return(Gamma)
  })

load_dt_pi <- nimbleFunction(
  run = function(xi=double(), nd=integer()) {
    returnType(double(1))
    pi <- rep(0,3+nd)
    pi[1] <- 1-xi
    pi[2] <- xi
    return(pi)
  })

load_dt_Pmats <- nimbleFunction(
  run = function(p=double(), nd=integer()) {
    returnType(double(2))
    Pmat <- nimMatrix(nrow=nd+3, ncol=2)
    Pmat[,1] <- 1
    Pmat[2,1] <- 1-p
    Pmat[2,2] <- p
    return(Pmat)
  })


