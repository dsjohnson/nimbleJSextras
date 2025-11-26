
#' @export
#' @import nimble
dnu <- nimble::nimbleFunction(
  run = function(x = double(2),
                 nUndet = integer(),
                 pi = double(1),        # Initial state probabilities
                 PArray = double(2),
                 GammaArray = double(3),   # Transition probabilities
                 len = integer(),
                 log = double(default=0)
  ) {
    returnType(double())
    out <- 0
    return(out)
  }
)

#' @export
#' @import nimble
rnu <- nimble::nimbleFunction(
  run = function(n = integer(default=1),
                 nUndet = integer(),
                 pi = double(1),        # Initial state probabilities
                 PArray = double(2),
                 GammaArray = double(3),   # Transition probabilities
                 len = integer(0)
  ) {
    # Number of states
    J <- dim(GammaArray)[1]
    K <- len

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=J)
    nu <- matrix(0,J,K)  # To store sampled states
    S <- matrix(0,J,J)

    # ---- Forward Pass ----
    Pdiag <- 1 - PArray[1,]
    alpha[1, ] <- init * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      Pdiag <- 1 - PArray[1,]
      alpha[t, ] <- (alpha[t-1, ] %*% GammaArray[, ,t-1]) * Pdiag
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    nu[,K] <- rmulti(1, nUndet, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      for(l in 1:J){
        trans_probs <- GammaArray[, l, rev_t] * alpha[rev_t, ]
        trans_probs <- trans_probs / sum(trans_probs)
        # Number of indiv in each state that transition to state "l"
        S[,l] <- rmulti(1, nu[l,rev_t+1], trans_probs)
      }
      # Sum across rows to get the number of indiv in each state at time K-t
      for(l in 1:J){
        nu[l,rev_t] <- sum(S[l,])
      }

    }
    returnType(double(2))
    return(nu)
  }
)
