#' Marked JS Observed History Conditional State Sampler Distributions
#'
#' @description
#' Custom NIMBLE distribution wrappers for \code{sample_det_ms}, \code{sample_det_binom},
#' \code{sample_det_pois}, and \code{sample_det_mark}. Provides dummy density functions
#' (returning constant log-density 0.0) and random simulation functions for conditional
#' latent state draws of observed individuals under various Jolly-Seber observation models.
#' @aliases dJSCond_ms dJSCond_binom dJSCond_pois dJSCond_mark rJSCond_ms rJSCond_binom rJSCond_pois rJSCond_mark
#' @name dJSCond

#' @param x Vector of sampled latent state histories of length K.
#' @param n Integer number of draws to simulate (defaults to 1).
#' @param obs Vector of observed capture history / counts / presence of length K.
#' @param init Vector of initial state distribution probabilities.
#' @param probObs 3D array of state-to-observation probabilities (n_states x O x K).
#' @param probTrans 3D array of state transition probability matrices ((J+2) x (J+2) x (K-1)).
#' @param prob 2D matrix of state-specific detection probabilities (K x J).
#' @param size Vector of binomial trial counts across occasions of length K.
#' @param rate 2D matrix of Poisson detection rates (K x J).
#' @param mark Vector of continuous mark observations of length K.
#' @param knots Vector of knot locations/means corresponding to each active state J.
#' @param sigma_knot Standard deviation around knot locations.
#' @param log Logical; if TRUE, returns log probability (always 0.0).
#'
#' @return \code{dJSCond_*} returns a scalar double (0.0).
#' \code{rJSCond_*} returns a numeric vector of sampled latent states of length K.
#'
#' -----------------------------------------------------------------------------
#' Multistate JS Conditional Sampler
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
dJSCond_ms <- nimbleFunction(
  run = function(x = double(1),
                 obs = double(1),
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 log = integer(default = 0)) {
    returnType(double(0))
    if (log) return(0.0)
    return(1.0)
  }
)

#' @export
#' @import nimble
#' @rdname dJSCond
rJSCond_ms <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 obs = double(1),
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3)) {
    returnType(double(1))
    sampled_z <- sample_det_ms(
      x = obs,
      init = init,
      probObs = probObs,
      probTrans = probTrans
    )
    return(sampled_z)
  }
)

#' -----------------------------------------------------------------------------
#' Binomial Observation JS Conditional Sampler (Robust Design)
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname dJSCond
dJSCond_binom <- nimbleFunction(
  run = function(x = double(1),
                 obs = double(1),
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 log = integer(default = 0)) {
    returnType(double(0))
    if (log) return(0.0)
    return(1.0)
  }
)

#' @export
#' @import nimble
#' @rdname dJSCond
rJSCond_binom <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 obs = double(1),
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3)) {
    returnType(double(1))
    sampled_z <- sample_det_binom(
      x = obs,
      init = init,
      prob = prob,
      size = size,
      probTrans = probTrans
    )
    return(sampled_z)
  }
)

#' -----------------------------------------------------------------------------
#' Poisson Observation JS Conditional Sampler (Continuous-time Robust Design)
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname dJSCond
dJSCond_pois <- nimbleFunction(
  run = function(x = double(1),
                 obs = double(1),
                 init = double(1),
                 rate = double(2),
                 probTrans = double(3),
                 log = integer(default = 0)) {
    returnType(double(0))
    if (log) return(0.0)
    return(1.0)
  }
)

#' @export
#' @import nimble
#' @rdname dJSCond
rJSCond_pois <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 obs = double(1),
                 init = double(1),
                 rate = double(2),
                 probTrans = double(3)) {
    returnType(double(1))
    sampled_z <- sample_det_pois(
      x = obs,
      init = init,
      rate = rate,
      probTrans = probTrans
    )
    return(sampled_z)
  }
)

#' -----------------------------------------------------------------------------
#' Bernoulli/Gaussian Marked Observation JS Conditional Sampler
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname dJSCond
dJSCond_mark <- nimbleFunction(
  run = function(x = double(1),
                 obs = double(1),
                 mark = double(1),
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 knots = double(1),
                 sigma_knot = double(),
                 log = integer(default = 0)) {
    returnType(double(0))
    if (log) return(0.0)
    return(1.0)
  }
)

#' @export
#' @import nimble
#' @rdname dJSCond
rJSCond_mark <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 obs = double(1),
                 mark = double(1),
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 knots = double(1),
                 sigma_knot = double()) {
    returnType(double(1))
    sampled_z <- sample_det_mark(
      x = obs,
      mark = mark,
      init = init,
      prob = prob,
      probTrans = probTrans,
      knots = knots,
      sigma_knot = sigma_knot
    )
    return(sampled_z)
  }
)
