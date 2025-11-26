
#' -----------------------------------------------------------------------------
#' Bernoulli (standard)
#' -----------------------------------------------------------------------------

#' @export
dJS <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 pstar = double(0),
                 piVector = double(1),##
                 PArray = double(2),
                 GammaArray = double(3),
                 len = integer(default=0),
                 log = integer(0, default = 0)) {
    len <- length(x)
    du <- dhmm_bern(x, piVector, PArray, GammaArray, 0, 1)
    logL <- du - log(pstar)
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
rJS <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 pstar = double(0, default=-1),
                 piVector = double(1),
                 PArray = double(2),
                 GammaArray = double(3),
                 len = integer()
  ) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rhmm_bern(n=1, piVector, PArray, GammaArray, len)
      ind <- ind+1
      if(all(x>0) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @export
calc_pstar <- nimble::nimbleFunction(
  run = function(piVector = double(1),
                 PArray = double(2),
                 GammaArray = double(3),
                 len = integer(0)) {
    zeros <- rep(0,len)
    pstar <- 1 - dhmm_bern(x=zeros, piVector, PArray, GammaArray, len, log=0)
    returnType(double(0))
    return(pstar)
  }
)

#' -----------------------------------------------------------------------------
#' Multistate
#' -----------------------------------------------------------------------------

#' @import nimble nimbleEcology
#' @export
dJS_ms <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 pstar = double(0),
                 piVector = double(1),##
                 PArray = double(3),
                 GammaArray = double(3),
                 len = integer(),
                 log = integer(0, default = 0)) {
    du <- dDHMMo(x, piVector, PArray, GammaArray, len, 0, 1)
    logL <- du - log(pstar)
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble nimbleEcology
#' @export
rJS_ms <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 pstar = double(0, default=-1),
                 piVector = double(1),##
                 PArray = double(3),
                 GammaArray = double(3),
                 len = integer()
  ) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rDHMMo(n=1, piVector, PArray, GammaArray, len, checkRowSums=0)
      ind <- ind+1
      if(!all(x==1) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @import nimble nimbleEcology
#' @export
calc_pstar_ms <- nimble::nimbleFunction(
  run = function(piVector = double(1),
                 PArray = double(3),
                 GammaArray = double(3),
                 len = integer()
                 ) {
    ones <- rep(1,len)
    pstar <- 1 - dDHMMo(ones, piVector, PArray, GammaArray, len, 0, 0)
    returnType(double(0))
    return(pstar)
  }
)


#' -----------------------------------------------------------------------------
#' Binomial robust design
#' -----------------------------------------------------------------------------

#' @export
dJS_rd <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 pstar = double(0),
                 piVector= double(1),##
                 PArray = double(2),
                 nSubOcc = double(1),
                 GammaArray = double(3),
                 len = integer(default=0),
                 log = integer(0, default = 0)) {
    len <- length(x)
    du <- dhmm_binom(x, piVector, PArray, nSubOcc, GammaArray, 0, log=1)
    logL <- du - log(pstar)
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
rJS_rd <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 pstar = double(0, default=-1),
                 piVector = double(1),
                 PArray = double(2),
                 nSubOcc = double(1),
                 GammaArray = double(3),
                 len = integer()) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rhmm_binom(n=1, piVector, PArray, nSubOcc, GammaArray, len)
      ind <- ind+1
      if(all(x>0) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @export
calc_pstar_rd <- nimble::nimbleFunction(
  run = function(piVector = double(1),
                 PArray = double(2),
                 nSubOcc = double(1),
                 GammaArray = double(3),
                 len = integer()) {
    zeros <- rep(0,len)
    pstar <- 1 - dhmm_binom(zeros, piVector, PArray, nSubOcc, GammaArray, len, 0)
    returnType(double(0))
    return(pstar)
  }
)
