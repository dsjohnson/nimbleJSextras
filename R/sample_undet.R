#' Sample states for undetected (not captured) individuals
#' @aliases sample_undet_binom sample_undet_ms sample_undet_pois
#' @name sample_undet
#' @param n Integer number of individuals to generate samples. This returns a
#' `n` x (number of occasions) matrix of sampled states for individuals that
#' were never captured.
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
#' @rdname sample_undet
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
#' @rdname sample_undet
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

