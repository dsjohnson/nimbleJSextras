#' Sample states for undetected (not captured) individuals
#' @aliases sample_det_binom sample_det_ms sample_det_pois
#' @name sample_det
#' @param x Observation sequence (capture-history) for an individual
#' @param init A vector of initial state probabilities
#' @param probObs Time-dependent 3D array of observation probabilities.  First
#' two dimensions of `probObs` are of size x (number of possible system states)
#' x (number of possible observation classes) and the third dimension is the
#' number of occasions. `probObs[i, j ,t]` is the probability that an individual
#' in the ith latent state is recorded as being in the jth detection state at
#' time t.
#' @param probTrans Time-dependent array of system state transition
#' probabilities. Dimension of `probTrans` is (number of possible system states)
#' x (number of possible system states) x (number of occasions).
#' `probTrans[i,j,t]` is the probability that an individual truly in state i at
#' time t will be in state j at time t+1.
#' @param prob For `dJS_binom` this is a 2d array that is (number of occasions)
#' x (number of states) which contains the probability of detection for each
#' state at time t. For a standard robust design there are 3 states: prior to
#' entry, entered, and exited. In this case the first and third column would
#' contain all 0s.
#' @param size A vector of the number of secondary occasions for each primary
#' occasion in a robust design.
#' @param rate A vector of detection rates for continuous-time secondary
#' occasions in `dJS_pois`. This should be proportional to the interval of time
#' each secondary occasion was executed.
#' @export
#' @import nimble
sample_det_ms <- nimble::nimbleFunction(
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
#' @import nimble
#' @rdname sample_det
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
#' @import nimble
#' @rdname sample_det
sample_det_pois <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 rate = double(2),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Number of states
    N <- dim(probTrans)[1]
    K <- length(x)  # Ensure len matches obs length

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    sampled_states <- numeric(K)  # To store sampled states

    # ---- Forward Pass ----
    Pdiag <- dpois(x[1], rate[1,])
    alpha[1, ] <- init * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
    for (t in 2:K) {
      Pdiag <- dpois(x[t], rate[t,])
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

