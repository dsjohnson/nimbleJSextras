

#' @import nimble nimbleEcology
#' @export
dJS_cat <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(3),
                 len = integer(),
                 pstar = double(0, default=-1),
                 checkRowSums = integer(0, default = 1),
                 log = integer(0, default = 0)) {
    du <- dDHMMo(x=x, probOb=probObs, probTrans=probTrans, init=init, len=len, log=1)
    if(pstar<0 | pstar>1){
      ones <- rep(1,len)
      dc <- 1 - dDHMMo(x=ones, probOb=probObs, probTrans=probTrans, init=init, len=len, log=0)
    } else{
      dc <- pstar
    }
    logL <- du - log(dc)
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble nimbleEcology
#' @export
rJS_cat <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(3),
                 len = integer(),
                 pstar = double(0, default=-1),
                 checkRowSums = integer(0, default = 1)
  ) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rDHMMo(n=1, init, probObs, probTrans, len, checkRowSums)
      ind <- ind+1
      if(!all(x==1) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @import nimble nimbleEcology
#' @export
pstar_cat <- nimble::nimbleFunction(
  run = function(init = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 len = integer(),
                 checkRowSums = integer(0, default = 1)) {
    ones <- rep(1,len)
    pstar <- 1 - dDHMMo(x=ones, probOb=probObs, probTrans=probTrans, init=init, len=len, log=0)
    returnType(double(0))
    return(pstar)
  }
)

#' -----------------------------------------------------------------------------
#' Bernoulli (standard)
#' -----------------------------------------------------------------------------

#' @export
dJS <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 prob = double(2),
                 probTrans = double(3),
                 pstar = double(0, default=-1),
                 len = integer(default=0),
                 log = integer(0, default = 0)) {
    len <- length(x)
    du <- dDHMMo_bern(x=x, init=init, prob=prob, probTrans=probTrans, len=len, log=1)
    if(pstar<0 | pstar>1){
      zeros <- rep(1,len)
      dc <- 1 - dDHMMo_bern(x=zeros, init=init, prob=prob, probTrans=probTrans, len=len, log=0)
    } else{
      dc <- pstar
    }
    logL <- du - log(dc)
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
rJS <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 pstar = double(0, default=-1),
                 len = integer(),
                 log = integer(0, default = 0)) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rDHMMo_bern(n=1, init, prob, probTrans, len, log)
      ind <- ind+1
      if(all(x>0) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @export
pstar_JS <- nimble::nimbleFunction(
  run = function(init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 len = integer(0, default=0)) {
    zeros <- rep(0,len)
    pstar <- 1 - dDHMMo_bern(x=zeros, init=init, prob=prob, probTrans=probTrans, len=len, log=0)
    returnType(double(0))
    return(pstar)
  }
)


#' -----------------------------------------------------------------------------
#' Binomial addition
#' -----------------------------------------------------------------------------

#' @export
dJS_binom <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 pstar = double(0, default=-1),
                 len = integer(default=0),
                 log = integer(0, default = 0)) {
    len <- length(x)
    du <- dDHMMo_binom(x=x, init=init, prob=prob, size=size, probTrans=probTrans, len=len, log=1)
    if(pstar<0 | pstar>1){
      zeros <- rep(1,len)
      dc <- 1 - dDHMMo_binom(x=zeros, init=init, prob=prob, size=size, probTrans=probTrans, len=len, log=0)
    } else{
      dc <- pstar
    }
    logL <- du - log(dc)
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
rJS_binom <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 pstar = double(0, default=-1),
                 len = integer(),
                 log = integer(0, default = 0)) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rDHMMo_binom(n=1, init, prob, size, probTrans, len, log)
      ind <- ind+1
      if(all(x>0) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @export
pstar_binom <- nimble::nimbleFunction(
  run = function(init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 len = integer(0, default=0)) {
    zeros <- rep(0,len)
    pstar <- 1 - dDHMMo_binom(x=zeros, init=init, prob=prob, size=size, probTrans=probTrans, len=len, log=0)
    returnType(double(0))
    return(pstar)
  }
)
