#' Convenience functions for loading parameter matrices for Jolly-Seber models
#'
#' @aliases load_Gamma load_Gamma_ms load_probObs load_probObs_ms
#' @name load_foo
#'
#' @param xi_tilde A double representing the recruitment/entry probability.
#' @param phi For \code{load_Gamma}: a scalar survival probability. For \code{load_Gamma_ms}: a numeric vector containing state-specific survival probabilities.
#' @param alpha A numeric vector specifying the initial distribution across alive states upon entry.
#' @param Psi A matrix of dimension S x S detailing the true movement transition probabilities.
#' @param p For \code{load_probObs}: a scalar capture probability. For \code{load_probObs_ms}: a numeric vector containing state-specific capture probabilities.
#'
#' @return A matrix representing the full state transition structure (for \code{load_Gamma} / \code{load_Gamma_ms})
#' or the full observation probability structure (for \code{load_probObs} / \code{load_probObs_ms}).
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


# #'@export
# #'@import nimble
# #'@rdname load_foo
# load_mix_Gamma <- nimbleFunction(
#   run = function(xi_tilde = double(1), phi = double(1)) {
#     returnType(double(2))
#
#     nmix <- length(phi)
#     idx <- c(1:(3*nmix))[rep(1:3, nmix)==1]-1
#     Gamma <- nimMatrix(0,3*nmix,3*nmix)
#     for(i in 1:nmix){
#       Gamma[idx[i]+1,idx[i]+1] <- 1 - xi_tilde[i]
#       Gamma[idx[i]+1,idx[i]+2] <- xi_tilde[i]
#       Gamma[idx[i]+2,idx[i]+2] <- phi[i]
#       Gamma[idx[i]+2,idx[i]+3] <- 1- phi[i]
#       Gamma[idx[i]+3,idx[i]+3] <- 1
#     }
#     return(Gamma)
#   })


# #'@export
# #'@import nimble
# #'@rdname load_foo
# load_mix_P <- nimbleFunction(
#   run = function(p = double(1)) {
#     returnType(double(2))
#
#     nmix <- length(p)
#     idx <- c(1:(3*nmix))[rep(1:3, nmix)==1]-1
#     P <- nimMatrix(0,3*nmix,2)
#     P[,1] <- 1
#     for(i in 1:nmix){
#       P[idx[i]+2,1] <- 1 - p[i]
#       P[idx[i]+2,2] <- p[i]
#     }
#     return(P)
#   })
