#' @export
#' @import nimble nimbleEcology
sample_undet_ms <- nimble::nimbleFunction(
  run = function(n = integer(0, default=1),
                 init = double(1),        # Initial state probabilities
                 probObs = double(3),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    J <- dim(probTrans)[1]
    K <- dim(probObs)[3]

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=J)
    sampled_states <- matrix(0,J,K)  # To store sampled states
    S <- matrix(0,J,J)

    # ---- Forward Pass ----
    alpha[1, ] <- init * probObs[, 1, 1]
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * probObs[, 1, t]
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    sampled_states[,K] <- rmulti(1, n, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      for(l in 1:J){
        trans_probs <- probTrans[, l, rev_t] * alpha[rev_t, ]
        trans_probs <- trans_probs / sum(trans_probs)
        # Number of indiv in each state that transition to state "l"
        S[,l] <- rmulti(1, sampled_states[l,rev_t+1], trans_probs)
      }
      # Sum across rows to get the number of indiv in each state at time K-t
      for(l in 1:J){
        sampled_states[l,rev_t] <- sum(S[l,])
      }

    }
    returnType(double(2))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' @export
#' @import nimble nimbleEcology
sample_undet_binom <- nimble::nimbleFunction(
  run = function(n = integer(0, default=1),
                 init = double(1),        # Initial state probabilities
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    J <- dim(probTrans)[1]
    K <- dim(prob)[1]

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=J)
    sampled_states <- matrix(0,J,K)  # To store sampled states
    S <- matrix(0,J,J)

    # ---- Forward Pass ----
    Pdiag <- dbinom(0, size[1], prob[1,])
    alpha[1, ] <- init * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      Pdiag <- dbinom(0, size[t], prob[t,])
      alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    sampled_states[,K] <- rmulti(1, n, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      for(l in 1:J){
        trans_probs <- probTrans[, l, rev_t] * alpha[rev_t, ]
        trans_probs <- trans_probs / sum(trans_probs)
        # Number of indiv in each state that transition to state "l"
        S[,l] <- rmulti(1, sampled_states[l,rev_t+1], trans_probs)
      }
      # Sum across rows to get the number of indiv in each state at time K-t
      for(l in 1:J){
        sampled_states[l,rev_t] <- sum(S[l,])
      }

    }
    returnType(double(2))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' @export
#' @import nimble nimbleEcology
sample_undet_pois <- nimble::nimbleFunction(
  run = function(n = integer(0, default=1),
                 init = double(1),        # Initial state probabilities
                 rate = double(2),
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    J <- dim(probTrans)[1]
    K <- dim(rate)[1]

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=J)
    sampled_states <- matrix(0,J,K)  # To store sampled states
    S <- matrix(0,J,J)

    # ---- Forward Pass ----
    Pdiag <- dpois(0, rate[1,])
    alpha[1, ] <- init * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
    for (t in 2:K) {
      Pdiag <- dpois(0, rate[t,])
      alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    sampled_states[,K] <- rmulti(1, n, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      for(l in 1:J){
        trans_probs <- probTrans[, l, rev_t] * alpha[rev_t, ]
        trans_probs <- trans_probs / sum(trans_probs)
        # Number of indiv in each state that transition to state "l"
        S[,l] <- rmulti(1, sampled_states[l,rev_t+1], trans_probs)
      }
      # Sum across rows to get the number of indiv in each state at time K-t
      for(l in 1:J){
        sampled_states[l,rev_t] <- sum(S[l,])
      }

    }
    returnType(double(2))  # Return sampled state sequence
    return(sampled_states)
  }
)

