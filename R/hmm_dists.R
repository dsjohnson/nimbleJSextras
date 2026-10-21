#' Hidden Markov Model Distributions for Binomial, Poisson, and Marked Capture Histories
#'
#' Density and random generation functions for Hidden Markov Model (HMM) distributions
#' with binomial, Poisson, or marked continuous observation distributions.
#'
#' @aliases dhmm_binom dhmm_pois rhmm_binom rhmm_pois dhmm_mark rhmm_mark
#' @name dhmm_js
#'
#' @param x Vector of observed counts or capture status for an individual across time (length `K`).
#' @param n Number of random draws, each returning a vector of length `len`.
#'   Currently only `n = 1` is supported.
#' @param init Vector of initial state probabilities at time `t = 1` of length `J+2`
#'   (where `J` is the number of intermediate active states, plus 1 initial unentered and 1 final exited state).
#' @param prob 2D array (`K x J`) containing detection probabilities for each
#'   intermediate/active state `J` at time `t`.
#' @param rate 2D array (`K x J`) containing Poisson detection rates for each
#'   intermediate/active state `J` at time `t`.
#' @param mark Vector of continuous observations/marks across time (length `K`).
#' @param knots Vector of knot locations/means corresponding to each intermediate state `J` (length `J`).
#' @param sigma_knot Standard deviation for the normal distribution around the knot locations.
#' @param size Vector of sizes (number of trials) for the binomial observation at each time (length `K`).
#' @param probTrans Time-dependent 3D array of system state transition probabilities
#'   with dimensions `(J+2) x (J+2) x K` (where `J+2` is the total number of system states and `K` is number of occasions).
#'   `probTrans[r, s, t]` is the probability of transitioning from state `r` at time `t` to state `s` at time `t+1`.
#' @param len Length of sequence `x` (`K`).
#' @param log Logical; if `TRUE`, log-probabilities/densities are returned.
#'
#' @return
#' For `dhmm_*`, the likelihood (or log-likelihood if `log = TRUE`) evaluated as a numeric scalar.
#' For `rhmm_*`, a numeric vector of generated values of length `len`.
#' Note that `rhmm_mark` currently returns a zero vector placeholder.
#'
#' -----------------------------------------------------------------------------
#' Binomial emission HMM distribution
#' -----------------------------------------------------------------------------

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
    returnType(double())
    pi <- init # State probabilities at time t=1
    logL <- 0
    K <- length(x)
    J <- dim(prob)[2]
    n_states <- J+2
    alpha <- pi
    pt <- nimNumeric(value = 0.0, length = n_states)
    for (t in 1:K) {
      ### Observation distribution -----------
      pt[2:(J+1)] <-  prob[t,]
      Pdiag <- dbinom(x[t], size[t], pt, log=0)
      ### ------------------------------------
      Zpi <- Pdiag * alpha
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t < K) alpha <- ((Zpi %*% probTrans[,,t]) / sumZpi)[1, ]
    }
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
    returnType(double(1))
    x <- numeric(len)
    n_states <- length(init)
    state <- rcat(1,init)
    for(t in 1:len){
      if(state==1 | state==n_states){
        p <- 0
      } else{
        p <- prob[t,state-1]
      }
      x[t] <- rbinom(1,size[t],p)
      if(t<len) state <- rcat(1,probTrans[state,,t])
    }
    return(x)
  }
)

#' -----------------------------------------------------------------------------
#' Poisson emission HMM distribution
#' -----------------------------------------------------------------------------

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
    returnType(double())
    pi <- init # State probabilities at time t=1
    logL <- 0
    K <- length(x)
    J <- dim(rate)[2]
    n_states <- J+2
    alpha <- pi
    lambda <- nimNumeric(value = 0.0, length = n_states)
    for (t in 1:K) {
      ### Observation distribution -----------
      lambda[2:(J+1)] <- rate[t,]
      Pdiag <- dpois(x[t], lambda, log=0)
      ### ------------------------------------
      Zpi <- Pdiag * alpha
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t < K) alpha <- ((Zpi %*% probTrans[,,t]) / sumZpi)[1, ]
    }
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
    returnType(double(1))
    x <- numeric(len)
    n_states <- length(init)
    state <- rcat(1,init)
    for(t in 1:len){
      if(state==1 | state==n_states){
        lambda <- 0
      } else{
        lambda <- rate[t,state-1]
      }
      x[t] <- rpois(1,lambda)
      if(t<len) state <- rcat(1,probTrans[state,,t])
    }
    return(x)
  }
)

#' -----------------------------------------------------------------------------
#' Marked (Bernoulli/Gaussian) emission HMM distribution
#' -----------------------------------------------------------------------------

#' @import nimble
#' @importFrom stats dnorm
#' @export
#' @rdname dhmm_js
dhmm_mark <- nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 mark = double(1), # mark vector -- can extend to 2d later
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 knots = double(1),
                 sigma_knot = double(),
                 len = integer(default=0),## length of x (needed as a separate param for rDHMMo_binom)
                 log = integer(default = 0)) {
    returnType(double())
    J <- length(knots)
    n_states <- J + 2
    pi <- init # State probabilities at time t=1
    logL <- 0
    len <- length(x)
    alpha <- pi
    dm <- nimNumeric(value = 1.0, length = n_states)
    pt <- nimNumeric(value = 0.0, length = n_states)
    for (t in 1:len) {
      ### Observation distribution -----------
      if(!is.na(mark[t])){
        dm[2:(J+1)] <- dnorm(mark[t], mean=knots, sd=sigma_knot, log=FALSE)
      } else {
        dm[2:(J+1)] <- 1.0
      }
      pt[2:(J+1)] <- prob[t,]
      Pdiag <- (1-x[t])*(1-pt) + x[t]*pt*dm
      ### ------------------------------------
      Zpi <- Pdiag * alpha
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t < len) alpha <- ((Zpi %*% probTrans[,,t]) / sumZpi)[1, ]
    }
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @import nimble
#' @export
#' @rdname dhmm_js
rhmm_mark <- nimbleFunction(
  run = function(n = integer(),
                 mark = double(1),
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 len = integer(0, default=0)
  ) {
    returnType(double(1))
    x <- numeric(len)
    return(x)
  }
)
