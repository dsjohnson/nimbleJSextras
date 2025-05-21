
#' @export
#' @import nimble nimbleEcology
smooth_HMM <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 probObs = double(2),     # Emission probabilities
                 probTrans = double(2)   # Transition probabilities
  ) {
    # Number of states
    N <- dim(probTrans)[1]
    K <- length(x)  # Ensure len matches obs length

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    beta <- matrix(0, nrow=K, ncol=N)
    gamma <- matrix(0, nrow=K, ncol=N)

    # ---- Forward Pass ----
    alpha[1, ] <- init * probObs[, x[1]]
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      alpha[t, ] <- (alpha[t-1, ] %*% probTrans) * probObs[, x[t]]
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Pass ----
    beta[K, ] <- alpha[K,]
    gamma[K, ] <- alpha[K, ]

    for (t in 1:(K-1)) {
      rev_t <- K-t
      beta[rev_t, ] <- probTrans %*% (probObs[, x[rev_t+1]] * beta[rev_t+1, ])
      beta[rev_t, ] <- beta[rev_t, ]/sum(beta[rev_t, ])
      gamma[rev_t, ] <- alpha[rev_t, ] * beta[rev_t, ]
      gamma[rev_t, ] <- gamma[rev_t, ]/sum(gamma[rev_t, ])
    }

    returnType(double(2))  # 2D matrix output
    return(gamma)
  }
)


#' @export
#' @import nimble nimbleEcology
smooth_HMMo <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 probObs = double(3),     # Emission probabilities
                 probTrans = double(2)   # Transition probabilities
  ) {
    # Number of states
    N <- dim(probTrans)[1]
    K <- length(x)  # Ensure len matches obs length

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    beta <- matrix(0, nrow=K, ncol=N)
    gamma <- matrix(0, nrow=K, ncol=N)

    # ---- Forward Pass ----
    alpha[1, ] <- init * probObs[, x[1], 1]
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      alpha[t, ] <- (alpha[t-1, ] %*% probTrans) * probObs[, x[t], t]
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Pass ----
    beta[K, ] <- alpha[K,]
    gamma[K, ] <- alpha[K, ]

    for (t in 1:(K-1)) {
      rev_t <- K-t
      beta[rev_t, ] <- probTrans %*% (probObs[, x[rev_t+1]] * beta[rev_t+1, ])
      beta[rev_t, ] <- beta[rev_t, ]/sum(beta[rev_t, ])
      gamma[rev_t, ] <- alpha[rev_t, ] * beta[rev_t, ]
      gamma[rev_t, ] <- gamma[rev_t, ]/sum(gamma[rev_t, ])
    }

    returnType(double(2))  # 2D matrix output
    return(gamma)
  }
)
