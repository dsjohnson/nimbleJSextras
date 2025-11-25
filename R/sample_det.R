


#' @export
#' @import nimble nimbleEcology
sample_det_cat <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 probObs = double(3),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    N <- dim(probTrans)[1]
    K <- length(x)  # Ensure len matches obs length

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    beta <- numeric(N)  # Only store last backward probabilities
    sampled_states <- numeric(K)  # To store sampled states

    # ---- Forward Pass ----
    alpha[1, ] <- init * probObs[, x[1], 1]
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * probObs[, x[t], t]
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    beta <- rep(1, N)  # Initialize backward probabilities
    sampled_states[K] <- rcat(1, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      trans_probs <- probTrans[, sampled_states[rev_t+1], rev_t] * alpha[rev_t, ]
      trans_probs <- trans_probs / sum(trans_probs)

      # Sample state directly
      sampled_states[rev_t] <- rcat(1, trans_probs)
    }

    returnType(double(1))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' @export
#' @import nimble nimbleEcology
sample_det_binom <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 prob = double(2),     # Emission probabilities
                 size = double(1),
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    N <- dim(probTrans)[1]
    K <- length(x)  # Ensure len matches obs length

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    sampled_states <- numeric(K)  # To store sampled states

    # ---- Forward Pass ----
    Pdiag <- dbinom(x[1], size[1], prob[1,])
    alpha[1, ] <- init * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      Pdiag <- dbinom(x[t], size[t], prob[t,])
      alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    sampled_states[K] <- rcat(1, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      trans_probs <- probTrans[, sampled_states[rev_t+1], rev_t] * alpha[rev_t, ]
      trans_probs <- trans_probs / sum(trans_probs)

      # Sample state directly
      sampled_states[rev_t] <- rcat(1, trans_probs)
    }

    returnType(double(1))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' @export
#' @import nimble nimbleEcology
sample_det <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 prob = double(2),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    N <- dim(probTrans)[1]
    K <- length(x)  # Ensure len matches obs length

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    sampled_states <- numeric(K)  # To store sampled states

    # ---- Forward Pass ----
    Pdiag <- dbinom(x[1], 1, prob[1,])
    alpha[1, ] <- init * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      Pdiag <- dbinom(x[t], 1, prob[t,])
      alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    sampled_states[K] <- rcat(1, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      trans_probs <- probTrans[, sampled_states[rev_t+1], rev_t] * alpha[rev_t, ]
      trans_probs <- trans_probs / sum(trans_probs)

      # Sample state directly
      sampled_states[rev_t] <- rcat(1, trans_probs)
    }

    returnType(double(1))  # Return sampled state sequence
    return(sampled_states)
  }
)



