#' -----------------------------------------------------------------------------
#' Binomial HMM distribution
#' -----------------------------------------------------------------------------

#' @import nimble
#' @export
dhmm_bern <- nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 pi = double(1),
                 PArray = double(2),
                 GammaArray = double(3),
                 len = integer(0, default=0),## length of x (needed as a separate param for rhmm_binom)
                 log = integer(0, default = 0)) {
    logL <- 0
    len <- length(x)
    for (t in 1:len) {
      ### Observation distribution -----------
      # Bernoulli
      Pdiag <- (1-x[t])*(1-PArray[t,]) + x[t]*PArray[t,]
      ### ------------------------------------
      Zpi <- Pdiag * pi
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t != len) pi <- ((Zpi %*% GammaArray[,,t]) / sumZpi)[1, ]
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble
#' @export
rhmm_bern <- nimbleFunction(
  run = function(n = integer(),
                 pi = double(1),
                 PArray = double(2),
                 GammaArray = double(3),
                 len = integer(0)## length of x (needed as a separate param for rDHMM)
                 ) {
    x <- numeric(len)
    state <- rcat(1,pi)
    for(t in 1:len){
      x[t] <- rbinom(1,1,PArray[t,state])
      if(t!=len) state <- rcat(1,GammaArray[state,,t])
    }
    returnType(double(1))
    return(x)
  }
)


#' -----------------------------------------------------------------------------
#' Binomial HMM distribution
#' -----------------------------------------------------------------------------

#' @import nimble
#' @export
dhmm_binom <- nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 pi = double(1),
                 PArray = double(2),
                 size = double(1),
                 GammaArray = double(3),
                 len = integer(0, default=0),## length of x (needed as a separate param for rhmm_binom)
                 log = integer(0, default = 0)) {
    logL <- 0
    len <- length(x)
    for (t in 1:len) {
      ### Observation distribution -----------
      Pdiag <- dbinom(x[t], size[t], PArray[t,])
      ### ------------------------------------
      Zpi <- Pdiag * pi
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t != len) pi <- ((Zpi %*% GammaArray[,,t]) / sumZpi)[1, ]
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble
#' @export
rhmm_binom <- nimbleFunction(
  run = function(n = integer(),
                 pi = double(1),
                 PArray = double(2),
                 size = double(1),
                 GammaArray = double(3),
                 len = integer(0)) {
    x <- numeric(len)
    state <- rcat(1,pi)
    for(t in 1:len){
      x[t] <- rbinom(1,size[t],PArray[t,state])
      if(t!=len) state <- rcat(1,GammaArray[state,,t])
    }
    returnType(double(1))
    return(x)
  }
)
