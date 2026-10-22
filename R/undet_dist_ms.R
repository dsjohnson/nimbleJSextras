#' Marked JS Unobserved Population Conditional State Sampler Distribution
#'
#' @description
#' Custom NIMBLE distribution wrappers for \code{sample_undet_ms}. Provides a dummy
#' density function and a random simulation function for conditional latent state draws
#' of unobserved individuals across time occasions.
#'
#' @param x Matrix of sampled latent states of dimension n_states x K.
#' @param n Integer number of draws to simulate (defaults to 1).
#' @param nu Scalar number of unobserved individuals in the superpopulation.
#' @param init Vector of initial state distribution probabilities.
#' @param probObs 3D array of state-to-observation probabilities of dimension n_states x O x K.
#' @param probTrans 3D array of state transition probability matrices of dimension n_states x n_states x (K-1).
#' @param log Logical; if TRUE, returns log probability (always 0.0).
#'
#' @return \code{dJSUndet_ms} returns a scalar double (0.0).
#' \code{rJSUndet_ms} returns a matrix of sampled state counts of dimension n_states x K.
#'
#' @export
#' @import nimble
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
rJSUndet_ms <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 nu = double(0),
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3)) {
    returnType(double(2))

    sampled_nu_t <- sample_undet_ms(
      n = nu,
      init = init,
      probObs = probObs,
      probTrans = probTrans
    )

    return(sampled_nu_t)
  }
)
