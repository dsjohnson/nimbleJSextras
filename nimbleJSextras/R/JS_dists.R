#' Jolly-Seber distributions For nimble models
#'
#' \code{dJS_ms} provides multistate Jolly-Seber (JS) distributions, \code{dJS_binom} provides binomial JS
#' distributions for standard JS and robust design JS models. \code{dJS_pois} provides JS distributions
#' for Poisson robust design, i.e., continuous-time secondary occasions. \code{pstar_**} functions
#' calculate the overall probability of detection and are necessary for the JS model distributions.
#'
#' @aliases dJS_ms dJS_binom dJS_pois dJS_mark rJS_ms rJS_binom rJS_pois rJS_mark pstar_ms pstar_binom pstar_pois pstar_mark
#' @name dJS
#'
#' @param x Vector of individual capture history (length `K`).
#' @param init Vector of initial state probabilities at time `t = 1` of length `J+2`
#'   (where `J` is the number of intermediate active states, plus 1 initial unentered and 1 final exited state).
#' @param probObs Time-dependent 3D array of observation probabilities. First
#' two dimensions of `probObs` are of size (number of possible system states)
#' `x` (number of possible observation classes) and the third dimension is the
#' number of occasions. `probObs[i, j, t]` is the probability that an individual
#' in the `i`th latent state is recorded as being in the `j`th detection state at
#' time `t`.
#' @param probTrans Time-dependent 3D array of system state transition probabilities
#'   with dimensions `(J+2) x (J+2) x K` (where `J+2` is the total number of system states and `K` is number of occasions).
#'   `probTrans[r, s, t]` is the probability of transitioning from state `r` at time `t` to state `s` at time `t+1`.
#' @param pstar Overall probability of detection. Calculated using the `pstar_**`
#' functions in nimble model code.
#' @param weight The number of times the capture history occurs in the full data
#' set of capture histories. Default is 1.
#' @param len Length of sequence `x` (`K`).
#' @param log Logical; if `TRUE`, log-probabilities/densities are returned.
#' @param n Number of random draws, each returning a vector of length `len`.
#' Currently only `n = 1` is supported, but the argument exists for
#' standardization of `rJS_**` functions.
#' @param prob For `dJS_binom`, a 2D array (`K x J`) containing detection probabilities
#' for each intermediate/active state `J` at time `t`.
#' @param size Vector of the number of secondary occasions for each primary
#' occasion in a robust design (length `K`).
#' @param rate 2D array (`K x J`) containing Poisson detection rates for each
#' intermediate/active state `J` at time `t`.
#' @param mark Vector of continuous observations/marks across time (length `K`).
#' @param knots Vector of knot locations/means corresponding to each intermediate state `J` (length `J`).
#' @param sigma_knot Standard deviation around knot locations.
#'
#' @return
#' For `dJS_*`, the likelihood (or log-likelihood if `log = TRUE`) evaluated as a numeric scalar.
#' For `rJS_*`, a numeric vector of generated capture histories of length `len`.
#' For `pstar_*`, a numeric scalar representing the overall detection probability.
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

# #' @import nimble nimbleEcology
# dJS_ms_n2ll <- nimble::nimbleFunction(
#   run = function(x = double(1),    ## Observed capture (state) history
#                  init = double(1),##
#                  probObs = double(3),
#                  probTrans = double(3),
#                  pstar = double(),
#                  weight = double(default=1),
#                  len = integer(default=1),
#                  log = integer(default = 0)
#   ) {
#     K <- length(x)
#     du <- dDHMMo(x, init, probObs, probTrans, K, 0, 1)
#     logn <- dpois(n, lambda*(1-pstar), log=1)/n
#     n2ll <- -2*weight*(du - log(pstar) + logn)
#     returnType(double(0))
#     return(n2ll)
#   }
# )


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
#' Poisson observations
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

#' -----------------------------------------------------------------------------
#' Bernoulli/Gaussian mark addition
#' -----------------------------------------------------------------------------

#' @import nimble
#' @export
#' @rdname dJS
dJS_mark <- nimble::nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 mark = double(1), # mark vector -- can extend to 2d later
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 pstar = double(),
                 knots = double(1),
                 sigma_knot = double(),
                 len = integer(default=0),## length of x (needed as a separate param for rDHMMo_binom)
                 log = integer(default = 0)) {
    returnType(double(0))
    K <- length(x)
    du <- dhmm_mark(x, mark, init, prob, probTrans, knots, sigma_knot, K, 1)
    logL <- du - log(pstar)
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble
#' @export
#' @rdname dJS
rJS_mark <- nimble::nimbleFunction(
  run = function(n = integer(),
                 mark = double(1),
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 pstar = double(),
                 knots = double(1),
                 sigma_knot = double(),
                 len = integer(default=0)
  ) {
    returnType(double(1))
    x <- nimNumeric(0, len)
    return(x)
  }
)

#' @import nimble
#' @export
#' @rdname dJS
pstar_mark <- nimble::nimbleFunction(
  run = function(init = double(1),
                 mark = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 knots = double(1),
                 sigma_knot = double(),
                 len = integer()) {
    zeros <- rep(0, len)
    pstar <- 1 - dhmm_mark(zeros, mark, init, prob, probTrans, knots, sigma_knot, len, 0)
    returnType(double(0))
    return(pstar)
  }
)
