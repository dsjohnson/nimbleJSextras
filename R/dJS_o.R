
#' @import nimble nimbleEcology
#' @export
dJS <- nimble::nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 init = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = integer(),
                 pstar = double(0, default=-1),
                 checkRowSums = integer(0, default = 1),
                 log = integer(0, default = 0)) {
    du <- dHMM(x=x, probOb=probObs, probTrans=probTrans, init=init, len=len, log=1)
    if(pstar<0 | pstar>1){
      ones <- rep(1,len)
      dc <- 1 - dHMM(x=ones, probOb=probObs, probTrans=probTrans, init=init, len=len, log=0)
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
rJS <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(2),
                 probTrans = double(2),
                 len = integer(),
                 pstar = double(0, default=-1),
                 checkRowSums = integer(0, default = 1)
  ) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rHMM(n=1, init, probObs, probTrans, len, checkRowSums)
      ind <- ind+1
      if(!all(x==1) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @import nimble nimbleEcology
#' @export
pstar <- nimble::nimbleFunction(
  run = function(init = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = integer(),
                 checkRowSums = integer(0, default = 1)) {
    ones <- rep(1,len)
    pstar <- 1 - dHMM(x=ones, probOb=probObs, probTrans=probTrans, init=init, len=len, log=0)
    returnType(double(0))
    return(pstar)
  }
)


#' -----------------------------------------------------------------------------
#' -----------------------------------------------------------------------------


#' @import nimble nimbleEcology
#' @export
dJS_o <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(2),
                 len = integer(),
                 pstar = double(0, default=-1),
                 checkRowSums = integer(0, default = 1),
                 log = integer(0, default = 0)) {
    du <- dHMMo(x=x, probOb=probObs, probTrans=probTrans, init=init, len=len, log=1)
    if(pstar<0 | pstar>1){
      ones <- rep(1,len)
      dc <- 1 - dHMMo(x=ones, probOb=probObs, probTrans=probTrans, init=init, len=len, log=0)
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
rJS_o <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(2),
                 len = integer(),
                 pstar = double(0, default=-1),
                 checkRowSums = integer(0, default = 1)
  ) {
    ind <- 0
    b <- 0
    while(b==0){
      x <- rHMMo(n=1, init, probObs, probTrans, len, checkRowSums)
      ind <- ind+1
      if(!all(x==1) | ind>10000) b=1
    }
    returnType(double(1))
    return(x)
  }
)

#' @import nimble nimbleEcology
#' @export
pstar_o <- nimble::nimbleFunction(
  run = function(init = double(1),
                 probObs = double(3),
                 probTrans = double(2),
                 len = integer(),
                 checkRowSums = integer(0, default = 1)) {
    ones <- rep(1,len)
    pstar <- 1 - dHMMo(x=ones, probOb=probObs, probTrans=probTrans, init=init, len=len, log=0)
    returnType(double(0))
    return(pstar)
  }
)


