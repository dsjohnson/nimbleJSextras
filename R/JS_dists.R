

#' @import nimble nimbleEcology
#' @export
dJS_ms <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(3),
                 pstar = double(),
                 weight = double(default=1),
                 len = integer(default=1),
                 log = integer(default = 0)) {
    K <- length(x)
    du <- dDHMMo(x, init, probObs, probTrans, K, 0, 1)
    # if(pstar<0 | pstar>1){
    #   ones <- rep(1,len)
    #   dc <- 1 - dDHMMo(x=ones, probOb=probObs, probTrans=probTrans, init=init, len=len, log=0)
    # } else{
    #   dc <- pstar
    # }
    logL <- weight*(du - log(pstar))
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble nimbleEcology
#' @export
rJS_ms <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(3),
                 pstar = double(default=0),
                 weight = double(default=1),
                 len = integer()
  ) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rDHMMo(n=1, init, probObs, probTrans, len, 0)
      ind <- ind+1
      if(!all(x==1) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @import nimble nimbleEcology
#' @export
pstar_ms <- nimble::nimbleFunction(
  run = function(init = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 len = integer()
  ) {
    ones <- rep(1,len)
    pstar <- 1 - dDHMMo(ones, init, probObs, probTrans, len, 0)
    returnType(double(0))
    return(pstar)
  }
)

#' @import nimble nimbleEcology
#' @export
dJS_ms_n2ll <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 n = integer(),
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(3),
                 pstar = double(),
                 lambda = double(),
                 weight = double(default=1)
  ) {
    K <- length(x)
    du <- dDHMMo(x, init, probObs, probTrans, K, 0, 1)
    # if(pstar<0 | pstar>1){
    #   ones <- rep(1,len)
    #   dc <- 1 - dDHMMo(x=ones, probOb=probObs, probTrans=probTrans, init=init, len=len, log=0)
    # } else{
    #   dc <- pstar
    # }
    logn <- dpois(n, lambda*(1-pstar), log=1)/n
    n2ll <- -2*weight*(du - log(pstar) + logn)
    returnType(double(0))
    return(n2ll)
  }
)


#' -----------------------------------------------------------------------------
#' Binomial addition
#' -----------------------------------------------------------------------------

#' @import nimble
#' @export
dJS_binom <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 pstar = double(),
                 weight = double(default=1),
                 len = integer(default=0),
                 log = integer(default = 0)) {
    K <- length(x)
    du <- dhmm_binom(x, init, prob, size, probTrans, K, 1)
    logL <- weight*(du - log(pstar))
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble
#' @export
rJS_binom <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 pstar = double(default=0),
                 weight = double(default=1),
                 len = integer()
  ) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rhmm_binom(n=1, init, prob, size, probTrans, len)
      ind <- ind+1
      if(all(x>0) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @import nimble
#' @export
pstar_binom <- nimble::nimbleFunction(
  run = function(init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 len = integer()) {
    zeros <- rep(0,len)
    pstar <- 1 - dhmm_binom(zeros, init, prob, size, probTrans, len, 0)
    returnType(double(0))
    return(pstar)
  }
)

#' -----------------------------------------------------------------------------
#' Poission observations
#' -----------------------------------------------------------------------------

#' @import nimble
#' @export
dJS_pois <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 rate = double(2),
                 probTrans = double(3),
                 pstar = double(),
                 weight = double(default=1),
                 len = integer(default=0),
                 log = integer(default = 0)) {
    K <- length(x)
    du <- dhmm_pois(x, init, rate, probTrans, K, 1)
    logL <- weight*(du - log(pstar))
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble
#' @export
rJS_pois <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 rate = double(2),
                 probTrans = double(3),
                 pstar = double(default=0),
                 weight = double(default=1),
                 len = integer()
  ) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rhmm_pois(n=1, init, rate, probTrans, len)
      ind <- ind+1
      if(all(x>0) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @import nimble
#' @export
pstar_pois <- nimble::nimbleFunction(
  run = function(init = double(1),
                 rate = double(2),
                 probTrans = double(3),
                 len = integer()) {
    zeros <- rep(0,len)
    pstar <- 1 - dhmm_pois(zeros, init, rate, probTrans, len, 0)
    returnType(double(0))
    return(pstar)
  }
)

