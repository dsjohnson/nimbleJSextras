#' Marked JS Unobserved Abundance Latent State Distributions
#'
#' @description
#' Custom NIMBLE distribution wrappers for sampling unobserved latent state histories
#' (nu_t) across time occasions for undetected individuals in Jolly-Seber models.
#' Provides dummy density functions (returning constant log-density 0.0) and random
#' simulation functions for unobserved abundance states.
#'
#' @aliases dJSUndet dJSUndet_ms rJSUndet_ms dJSUndet_binom rJSUndet_binom dJSUndet_pois rJSUndet_pois dJSUndet_mark rJSUndet_mark
#' @name dJSUndet
#'
#' @param x Matrix of sampled unobserved latent states across intermediate states and occasions (dimension J x K).
#' @param n Integer number of draws to simulate (defaults to 1).
#' @param nu Total number of unobserved individuals in the superpopulation.
#' @param init Vector of initial state distribution probabilities.
#' @param probObs 3D array of state-to-observation probabilities (n_states x O x K).
#' @param probTrans 3D array of state transition probability matrices ((J+2) x (J+2) x (K-1)).
#' @param prob 2D matrix of state-specific detection probabilities (K x J).
#' @param size Vector of binomial trial counts across occasions of length K.
#' @param rate 2D matrix of Poisson detection rates (K x J).
#' @param log Logical; if TRUE, returns log probability (always 0.0).
#'
#' @return \code{dJSUndet_*} returns a scalar double (0.0).
#' \code{rJSUndet_*} returns a 2D numeric matrix of sampled unobserved latent state counts across occasions.
#'
#' -----------------------------------------------------------------------------
#' Multistate JS Unobserved Sampler
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname dJSUndet
dJSUndet_ms <- nimbleFunction(
  run = function(x = double(2),
                 nu = double(0),
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
#' @rdname dJSUndet
rJSUndet_ms <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 nu = double(0),
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3)) {
    returnType(double(2))
    sampled_nu <- sample_undet_ms(
      n = nu,
      init = init,
      probObs = probObs,
      probTrans = probTrans
    )
    return(sampled_nu)
  }
)

#' -----------------------------------------------------------------------------
#' Binomial Observation JS Unobserved Sampler (Robust Design)
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname dJSUndet
dJSUndet_binom <- nimbleFunction(
  run = function(x = double(2),
                 nu = double(0),
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
#' @rdname dJSUndet
rJSUndet_binom <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 nu = double(0),
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3)) {
    returnType(double(2))
    sampled_nu <- sample_undet_binom(
      n = nu,
      init = init,
      prob = prob,
      size = size,
      probTrans = probTrans
    )
    return(sampled_nu)
  }
)

#' -----------------------------------------------------------------------------
#' Poisson Observation JS Unobserved Sampler (Continuous-time Robust Design)
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname dJSUndet
dJSUndet_pois <- nimbleFunction(
  run = function(x = double(2),
                 nu = double(0),
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
#' @rdname dJSUndet
rJSUndet_pois <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 nu = double(0),
                 init = double(1),
                 rate = double(2),
                 probTrans = double(3)) {
    returnType(double(2))
    sampled_nu <- sample_undet_pois(
      n = nu,
      init = init,
      rate = rate,
      probTrans = probTrans
    )
    return(sampled_nu)
  }
)

#' -----------------------------------------------------------------------------
#' Bernoulli/Gaussian Marked Observation JS Unobserved Sampler
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname dJSUndet
dJSUndet_mark <- nimbleFunction(
  run = function(x = double(2),
                 nu = double(0),
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 log = integer(default = 0)) {
    returnType(double(0))
    if (log) return(0.0)
    return(1.0)
  }
)

#' @export
#' @import nimble
#' @rdname dJSUndet
rJSUndet_mark <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 nu = double(0),
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3)) {
    returnType(double(2))
    sampled_nu <- sample_undet_mark(
      n = nu,
      init = init,
      prob = prob,
      probTrans = probTrans
    )
    return(sampled_nu)
  }
)
