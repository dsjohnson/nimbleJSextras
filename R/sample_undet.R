#' Sample latent states for undetected (unobserved) individuals in Jolly-Seber models
#'
#' Forward-filtering backward-sampling functions to simulate state trajectories
#' for individuals that were never observed in the capture-recapture study.
#'
#' @aliases sample_undet sample_undet_ms sample_undet_binom sample_undet_pois sample_undet_mark
#' @name sample_undet
#'
#' @param n Integer number of undetected individuals to simulate.
#' @param init Vector of initial state probabilities at time `t = 1` of length `J+2`
#'   (where `J` is the number of intermediate active states, plus 1 initial unentered and 1 final exited state).
#' @param probObs Time-dependent 3D array of observation probabilities. First
#'   two dimensions of `probObs` are of size `(J+2) x O` (`O` = number of observation states)
#'   and the third dimension is `K`.
#' @param probTrans Time-dependent 3D array of system state transition probabilities
#'   with dimensions `(J+2) x (J+2) x K`.
#'   `probTrans[r, s, t]` is the probability of transitioning from state `r` at time `t` to state `s` at time `t+1`.
#' @param prob 2D array (`K x J`) containing detection probabilities for each
#'   intermediate/active state `J` at time `t`.
#' @param size Vector of sizes (number of trials) for the binomial observation at each primary occasion (length `K`).
#' @param rate 2D array (`K x J`) containing Poisson detection rates for each
#'   intermediate/active state `J` at time `t`.
#'
#' @return A matrix of dimension `(J+2) x K` containing the count of undetected individuals
#' in each system state across all occasions `K`.
#'
#' -----------------------------------------------------------------------------
#' Multistate JS
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble nimbleEcology
#' @rdname sample_undet
sample_undet_ms <- nimble::nimbleFunction(
  run = function(n = integer(0, default=1),
                 init = double(1),        # Initial state probabilities
                 probObs = double(3),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Dimensions
    n_states <- dim(probTrans)[1]
    K <- dim(probObs)[3]
    # Forward probabilities
    alpha <- matrix(0, nrow=K, ncol=n_states)
    sampled_states <- matrix(0,n_states,K)  # To store sampled states
    S <- matrix(0,n_states,n_states)
    # ---- Forward Pass ----
    for (t in 1:K) {
      if(t == 1){
        alpha[t, ] <- init * probObs[, 1, t]
      } else{
        alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * probObs[, 1, t]
      }
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }
    # ---- Backward Sampling ----
    sampled_states[,K] <- rmulti(1, n, alpha[K, ])  # Sample final state
    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      for(k in 1:n_states){
        back_probs <- probTrans[, k, rev_t] * alpha[rev_t, ]
        back_probs <- back_probs / sum(back_probs)
        # Number of indiv in each state that transition to state "l"
        S[,k] <- rmulti(1, sampled_states[k,rev_t+1], back_probs)
      }
      # Sum across rows to get the number of indiv in each state at time K-t
      for(k in 1:n_states){
        sampled_states[k,rev_t] <- sum(S[k,])
      }
    }
    returnType(double(2))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' -----------------------------------------------------------------------------
#' Binomial
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble nimbleEcology
#' @rdname sample_undet
sample_undet_binom <- nimble::nimbleFunction(
  run = function(n = integer(0, default=1),
                 init = double(1),        # Initial state probabilities
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    n_states <- length(init)
    K <- dim(prob)[1]
    J <- dim(prob)[2]
    # Forward probabilities
    alpha <- matrix(0, nrow=K, ncol=n_states)
    sampled_states <- matrix(0,n_states,K)  # To store sampled states
    S <- matrix(0,n_states,n_states)
    Pdiag <- nimNumeric(value = 1.0, length = n_states)
    # ---- Forward Pass ----
    for (t in 1:K) {
      Pdiag[2:(J+1)] <- (1-prob[t,])^size[t]
      if(t == 1){
        alpha[t, ] <- init * Pdiag
      } else{
        alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      }
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    sampled_states[,K] <- rmulti(1, n, alpha[K, ])  # Sample final state
    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      for(k in 1:n_states){
        back_probs <- probTrans[, k, rev_t] * alpha[rev_t, ]
        back_probs <- back_probs / sum(back_probs)
        # Number of indiv in each state that transition to state "l"
        S[,k] <- rmulti(1, sampled_states[k,rev_t+1], back_probs)
      }
      # Sum across rows to get the number of indiv in each state at time K-t
      for(k in 1:n_states){
        sampled_states[k,rev_t] <- sum(S[k,])
      }
    }
    returnType(double(2))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' -----------------------------------------------------------------------------
#' Poisson
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble nimbleEcology
#' @rdname sample_undet
sample_undet_pois <- nimble::nimbleFunction(
  run = function(n = integer(0, default=1),
                 init = double(1),        # Initial state probabilities
                 rate = double(2),
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    n_states <- length(init)
    K <- dim(rate)[1]
    J <- dim(rate)[2]
    # Forward probabilities
    alpha <- matrix(0, nrow=K, ncol=n_states)
    sampled_states <- matrix(0,n_states,K)  # To store sampled states
    S <- matrix(0,n_states,n_states)
    Pdiag <- nimNumeric(value = 1.0, length = n_states)
    # ---- Forward Pass ----
    for (t in 1:K) {
      Pdiag[2:(J+1)] <- exp(-rate[t,])
      if(t == 1){
        alpha[t, ] <- init * Pdiag
      } else {
        alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      }
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }
    # ---- Backward Sampling ----
    sampled_states[,K] <- rmulti(1, n, alpha[K, ])  # Sample final state
    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      for(k in 1:n_states){
        back_probs <- probTrans[, k, rev_t] * alpha[rev_t, ]
        back_probs <- back_probs / sum(back_probs)
        # Number of indiv in each state that transition to state "l"
        S[,k] <- rmulti(1, sampled_states[k,rev_t+1], back_probs)
      }
      # Sum across rows to get the number of indiv in each state at time K-t
      for(k in 1:n_states){
        sampled_states[k,rev_t] <- sum(S[k,])
      }
    }
    returnType(double(2))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' -----------------------------------------------------------------------------
#' Bernoulli/ Gaussian mark
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble nimbleEcology
#' @rdname sample_undet
sample_undet_mark <- nimble::nimbleFunction(
  run = function(n = integer(0, default=1),
                 init = double(1),        # Initial state probabilities
                 prob = double(2),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
  ) {
    n_states <- length(init)
    K <- dim(prob)[1]
    J <- dim(prob)[2]
    # Forward probabilities
    alpha <- matrix(0, nrow=K, ncol=n_states)
    sampled_states <- matrix(0,n_states,K)  # To store sampled states
    S <- matrix(0,n_states,n_states)
    Pdiag <- nimNumeric(value = 1.0, length = n_states)
    # ---- Forward Pass ----
    for (t in 1:K) {
      Pdiag[2:(J+1)] <- 1-prob[t,]
      if(t == 1){
        alpha[t, ] <- init * Pdiag
      } else{
        alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      }
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }
    # ---- Backward Sampling ----
    sampled_states[,K] <- rmulti(1, n, alpha[K, ])  # Sample final state
    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      for(k in 1:n_states){
        back_probs <- probTrans[, k, rev_t] * alpha[rev_t, ]
        back_probs <- back_probs / sum(back_probs)
        # Number of indiv in each state that transition to state "l"
        S[,k] <- rmulti(1, sampled_states[k,rev_t+1], back_probs)
      }
      # Sum across rows to get the number of indiv in each state at time K-t
      for(k in 1:n_states){
        sampled_states[k,rev_t] <- sum(S[k,])
      }
    }
    returnType(double(2))  # Return sampled state sequence
    return(sampled_states)
  }
)
