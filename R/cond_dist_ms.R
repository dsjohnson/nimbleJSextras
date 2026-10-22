#' Marked JS Observed History Conditional State Sampler Distribution
#'
#' @description
#' Custom NIMBLE distribution wrappers for \code{sample_det_ms}. Provides a dummy
#' density function and a random simulation function for conditional latent state draws
#' of observed individuals.
#'
#' @param x Vector of sampled latent state histories of length K.
#' @param n Integer number of draws to simulate (defaults to 1).
#' @param obs Vector of observed capture history of length K.
#' @param init Vector of initial state distribution probabilities.
#' @param probObs 3D array of state-to-observation probabilities of dimension n_states x O x K.
#' @param probTrans 3D array of state transition probability matrices of dimension n_states x n_states x (K-1).
#' @param log Logical; if TRUE, returns log probability (always 0.0).
#'
#' @return \code{dJSCond_ms} returns a scalar double (0.0).
#' \code{rJSCond_ms} returns a vector of sampled latent states of length K.
#'
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
