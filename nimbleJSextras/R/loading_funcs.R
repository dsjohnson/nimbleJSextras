#' Convenience functions for loading transition and observation parameter matrices for Jolly-Seber models
#'
#' Functions to generate the full transition matrix `Gamma` or observation probability
#' matrix `P` for single-state and multistate Jolly-Seber models.
#'
#' @aliases load_foo load_Gamma load_Gamma_ms load_probObs load_probObs_ms
#' @name load_foo
#'
#' @param xi_tilde Scalar representing recruitment/entry probability.
#' @param phi For \code{load_Gamma}: a scalar survival probability.
#'   For \code{load_Gamma_ms}: a numeric vector of length `J` containing state-specific survival probabilities.
#' @param alpha A numeric vector of length `J` specifying the initial distribution across active states upon entry.
#' @param Psi Matrix of dimension `J x J` detailing movement transition probabilities between active states.
#' @param p For \code{load_probObs}: a scalar capture probability.
#'   For \code{load_probObs_ms}: a numeric vector of length `J` containing state-specific capture probabilities.
#'
#' @return
#' For single-state models (\code{load_Gamma}, \code{load_probObs}):
#' \itemize{
#'   \item \code{load_Gamma}: A `3 x 3` state-transition matrix.
#'   \item \code{load_probObs}: A `3 x 2` observation probability matrix.
#' }
#' For multistate models (\code{load_Gamma_ms}, \code{load_probObs_ms}):
#' \itemize{
#'   \item \code{load_Gamma_ms}: A `(J+2) x (J+2)` state-transition matrix (where `J` is the number of active states).
#'   \item \code{load_probObs_ms}: A `(J+2) x 2` observation probability matrix.
#' }
#'
#' @import nimble
NULL

#' @rdname load_foo
#' @export
#' @import nimble
load_Gamma <- nimbleFunction(
  run = function(xi_tilde = double(), phi = double()) {
    returnType(double(2))
    Gamma <- nimMatrix(0,3,3)
    Gamma[1,1] <- 1 - xi_tilde
    Gamma[1,2] <- xi_tilde
    Gamma[2,2] <- phi
    Gamma[2,3] <- 1- phi
    Gamma[3,3] <- 1
    return(Gamma)
  })

#' @rdname load_foo
#' @export
#' @import nimble
load_probObs <- nimbleFunction(
  run = function(p = double()) {
    returnType(double(2))
    P <- nimMatrix(0,3,2)
    P[1,1] <- 1
    P[1,2] <- 0
    P[2,1] <- 1-p
    P[2,2] <- p
    P[3,1] <- 1
    P[3,2] <- 0
    return(P)
  })


#' @rdname load_foo
#' @export
#' @import nimble
load_Gamma_ms <- nimbleFunction(
  run = function(
    xi_tilde = double(), alpha=double(1), phi = double(1), Psi=double(2)
  ) {
    returnType(double(2))
    S <- dim(Psi)[1]
    Gamma <- nimMatrix(0,S+2,S+2)
    Gamma[1,1] <- 1 - xi_tilde
    for(j in 1:S){
      Gamma[1,j+1] <- xi_tilde*alpha[j]
      for(k in 1:S) Gamma[j+1,k+1] <- phi[j]*Psi[j,k]
      Gamma[j+1,S+2] <- 1- phi[j]
    }
    Gamma[S+2,S+2] <- 1
    return(Gamma)
  })

#' @rdname load_foo
#' @export
#' @import nimble
load_probObs_ms <- nimbleFunction(
  run = function(p = double(1)) {
    returnType(double(2))
    S <- length(p)
    P <- nimMatrix(0,S+2,2)
    P[1,1] <- 1
    P[1,2] <- 0
    for(k in 1:S){
      P[k+1,1] <- 1-p[k]
      P[k+1,2] <- p[k]
    }
    P[S+2,1] <- 1
    P[S+2,2] <- 0
    return(P)
  })
