#' Hidden Markov Model distribution for binomial and Poisson capture histories
#' @aliases dhmm_binom dhmm_pois rhmm_binom rhmm_pois
#' @name dhmm_js
#'
#' @param x Vector of observed counts for an individual
#' @param init Vector of initial state probabilities.
#' @param prob 2d array that is (length of x)
#' x (number of states) which contains the probability of detection for each
#' state at time t.
#' @param size Vector of sizes for the binomial observation at each time.
#' @param probTrans Time-dependent array of system state transition
#' probabilities. Dimension of `probTrans` is (number of possible system states)
#' x (number of possible system states) x (number of occasions).
#' `probTrans[i,j,t]` is the probability that an individual truly in state i at
#' time t will be in state j at time t+1.
#' @param len Length of `x`
#' @param log Logical. Return log density
#' @param n number of random draws, each returning a vector of length `len`.
#' Currently only `n = 1` is supported, but the argument exists for
#' standardization of `rhmm_**` functions.
#' @param rate A vector of rates for each time in the `dhmm_pois` distribution function.

#' @import nimble
#' @export
dhmm_binom <- nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 len = integer(default=0),## length of x (needed as a separate param for rDHMMo_binom)
                 log = integer(default = 0)) {
    pi <- init # State probabilities at time t=1
    logL <- 0
    len <- length(x)
    alpha <- pi
    for (t in 1:len) {
      ### Observation distribution -----------
      Pdiag <- dbinom(x[t], size[t], prob[t,], log=0)
      ### ------------------------------------
      Zpi <- Pdiag * alpha
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t < len) alpha <- ((Zpi %*% probTrans[,,t]) / sumZpi)[1, ]
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble
#' @export
#' @rdname dhmm_js
rhmm_binom <- nimbleFunction(
  run = function(n = integer(),
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 len = integer(0, default=0)
                 ) {
    x <- numeric(len)
    state <- rcat(1,init)
    for(t in 1:len){
      x[t] <- rbinom(1,size[t],prob[t,state])
      if(t<len) state <- rcat(1,probTrans[state,,t])
    }
    returnType(double(1))
    return(x)
  }
)

#' @import nimble
#' @export
#' @rdname dhmm_js
dhmm_pois <- nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 init = double(1),
                 rate = double(2),
                 probTrans = double(3),
                 len = integer(default=0),## length of x (needed as a separate param for rDHMMo_binom)
                 log = integer(default = 0)) {
    pi <- init # State probabilities at time t=1
    logL <- 0
    len <- length(x)
    alpha <- pi
    for (t in 1:len) {
      ### Observation distribution -----------
      Pdiag <- dpois(x[t], rate[t,], log=0)
      ### ------------------------------------
      Zpi <- Pdiag * alpha
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t < len) alpha <- ((Zpi %*% probTrans[,,t]) / sumZpi)[1, ]
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble
#' @export
#' @rdname dhmm_js
rhmm_pois <- nimbleFunction(
  run = function(n = integer(),
                 init = double(1),
                 rate = double(2),
                 probTrans = double(3),
                 len = integer(0, default=0)
  ) {
    x <- numeric(len)
    state <- rcat(1,init)
    for(t in 1:len){
      x[t] <- rpois(1,rate[t,state])
      if(t<len) state <- rcat(1,probTrans[state,,t])
    }
    returnType(double(1))
    return(x)
  }
)

