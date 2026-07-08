
#' Jolly-Seber distributions For nimble models
#'
#' \code{dJS_ms} provides multistate Jolly-Seber (JS) distributions, \code{dJS_binom} provides binomial JS
#' distributions for standard JS and robust design JS models. \code{dJS_pois} provides JS distributions
#' for Poisson robust design, i.e., continuous-time secondary occasions. \code{pstar_**} functions
#' calculate the overall probability of detection and are necessary for the JS model distributions.
#'
#' @aliases dJS_ms dJS_binom dJS_pois rJS_ms rJS_binom rJS_pois pstar_ms pstar_binom pstar_pois
#' @name dJS
#'
#' @param x In dividual capture history
#' @param init A vector of initial state probability vectors.
#' @param probObs Time-dependent 3D array of observation probabilities.  First
#' two dimensions of `probObs` are of size x (number of possible system states)
#' x (number of possible observation classes) and the third dimension is the
#' number of occasions. `probObs[i, j ,t]` is the probability that an individual
#' in the ith latent state is recorded as being in the jth detection state at
#' time t.
#' @param probTrans Time-dependent array of system state transition
#' probabilities. Dimension of `probTrans` is (number of possible system states)
#' x (number of possible system states) x (number of occasions).
#' `probTrans[i,j,t]` is the probability that an individual truly in state i at
#' time t will be in state j at time t+1.
#' @param pstar Overall probability of detection. Calculated using the `pstar_**`
#' functions in nimble model code.
#' @param weight The number of times the capture history occurs in the full data
#' set of capture histories.
#' @param len Length of `x`
#' @param log Logical. Return log-distribution.
#' @param n number of random draws, each returning a vector of length `len`.
#' Currently only `n = 1` is supported, but the argument exists for
#' standardization of `rJS_**` functions.
#' @param prob For `dJS_binom` this is a 2d array that is (number of occasions)
#' x (number of states) which contains the probability of detection for each
#' state at time t. For a standard robust design there are 3 states: prior to
#' entry, entered, and exited. In this case the first and third column would
#' contain all 0s.
#' @param size A vector of the number of secondary occasions for each primary
#' occasion in a robust design.
#' @param rate A vector of detection rates for continuous-time secondary
#' occasions in `dJS_pois`. This should be proportional to the interval of time
#' each secondary occasion was executed.
#'
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
    logL <- weight*(du - log(pstar))
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble nimbleEcology
#' @export
#' @rdname dJS
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
#' @rdname dJS
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
#' @rdname dJS
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
#' @rdname dJS
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
#' @rdname dJS
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
#' @rdname dJS
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
#' @rdname dJS
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
#' @rdname dJS
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

