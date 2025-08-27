#' @export
#' @import nimble nimbleEcology
sample_state_undet_Do <- nimble::nimbleFunction(
  run = function(n = integer(0),
                 init = double(1),        # Initial state probabilities
                 probObs = double(3),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
                 len = integer(0)
  ) {
    # Number of states
    J <- dim(probTrans)[1]
    K <- len

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=J)
    sampled_states <- matrix(0,J,K)  # To store sampled states

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

        # Sample state directly
        !!!!!!!
        sampled_states[j,rev_t] <- rmulti(1, sampled_states[l,rev_t+1], trans_probs)
      }
    }

    returnType(double(1))  # Return sampled state sequence
    return(sampled_states)
  }
)


#' @import nimble
#' @export
sample_state_undet_Do <- nimble::nimbleFunction(
  run = function(
    nu = integer(),           # Total number of undetected indiv.
    init = double(1),        # Initial state probabilities
    probTrans = double(3),    # Transition probabilities
    probObs = double(3),
    len = integer(0, default=0)
  ) {
    J <- dim(probTrans)[1]
    K <- dim(probTrans)[3] + 1
    nu_t <- matrix(0,nrow=J, ncol=K)
    Tu <- array(0,dim = c(J,J,K-1))
    nu_t[1:J,1] <- rmulti(1, nu, init[1:J]*probObs[1:J,1,1])
    for(t in 2:K){
      for(l in 1:J){
        Tu[l,1:J,t-1] <- rmulti(1,  nu_t[l,t-1], probTrans[l,1:J,t-1]*probObs[1:J,1,t])
        nu_t[l,t] <- sum(Tu[1:J,l,t-1])
      }
    }
    returnType(double(2))  # Return sampled state sequence
    return(nu_t)
  }
)
