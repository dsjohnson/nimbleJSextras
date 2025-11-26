
###
### Sample true states for observed individuals
###

## Standard JS model ###########################################################

#' @export
dState <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 PiVector = double(1),##
                 PArray = double(2),
                 GammaArray = double(3),
                 x_cond = double(1),
                 log = integer(0, default = 0)) {
    logL <- 0
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
rState <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 PiVector = double(1),
                 PArray = double(2),
                 GammaArray = double(3),
                 x_cond = double(1)
  ) {

    x <- x_cond
    N <- dim(GammaArray)[1]
    K <- length(x)

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    z <- numeric(K)  # To store sampled states

    # ---- Forward Pass ----
    Pdiag <- (1-x[1])*(1-PArray[1,]) + x[t]*PArray[1,]
    alpha[1, ] <- PiVector * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      Pdiag <- (1-x[t])*(1-PArray[t,]) + x[t]*PArray[t,]
      alpha[t, ] <- (alpha[t-1, ] %*% GammaArray[, ,t-1]) * Pdiag
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    z[K] <- rcat(1, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      trans_probs <- GammaArray[, z[rev_t+1], rev_t] * alpha[rev_t, ]
      trans_probs <- trans_probs / sum(trans_probs)

      # Sample state directly
      z[rev_t] <- rcat(1, trans_probs)
    }

    returnType(double(1))
    return(z)
  }
)

## Multistate model ###########################################################

#' @export
dState_ms <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 PiVector = double(1),##
                 PArray = double(2),
                 GammaArray = double(3),
                 x_cond = double(1),
                 log = integer(0, default = 0)) {
    logL <- 0
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
rState_ms <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 PiVector = double(1),
                 PArray = double(2),
                 GammaArray = double(3),
                 x_cond = double(1)
  ) {

    x <- x_cond
    N <- dim(GammaArray)[1]
    K <- length(x)  # Ensure len matches obs length

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    beta <- numeric(N)  # Only store last backward probabilities
    z <- numeric(K)  # To store sampled states

    # ---- Forward Pass ----
    alpha[1, ] <- init * PAarray[, x[1], 1]
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      alpha[t, ] <- (alpha[t-1, ] %*% GammaArray[, ,t-1]) * PArray[, x[t], t]
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    beta <- rep(1, N)  # Initialize backward probabilities
    z[K] <- rcat(1, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      trans_probs <- GammaArray[, z[rev_t+1], rev_t] * alpha[rev_t, ]
      trans_probs <- trans_probs / sum(trans_probs)

      # Sample state directly
      z[rev_t] <- rcat(1, trans_probs)
    }

    returnType(double(1))  # Return sampled state sequence
    return(z)
  }
)

## Robust design (Binomial) model ###########################################################

#' @export
dState_rd <- nimble::nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 PiVector = double(1),##
                 PArray = double(2),
                 nSubOcc = double(1),
                 GammaArray = double(3),
                 x_cond = double(1),
                 log = integer(0, default = 0)) {
    logL <- 0
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
rState_rd <- nimble::nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 PiVector = double(1),
                 PArray = double(2),
                 nSubOcc = double(1),
                 GammaArray = double(3),
                 x_cond = double(1)
  ) {

    x <- x_cond
    N <- dim(GammaArray)[1]
    K <- length(x)

    # Forward (α) and Backward (β) matrices
    alpha <- matrix(0, nrow=K, ncol=N)
    z <- numeric(K)  # To store sampled states

    # ---- Forward Pass ----
    Pdiag <- dbinom(x[1], nSubOcc[1], PArray[1,])
    alpha[1, ] <- PiVector * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

    for (t in 2:K) {
      Pdiag <- dbinom(x[t], nSubOcc[t], PArray[t,])
      alpha[t, ] <- (alpha[t-1, ] %*% GammaArray[, ,t-1]) * Pdiag
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }

    # ---- Backward Sampling ----
    z[K] <- rcat(1, alpha[K, ])  # Sample final state

    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      trans_probs <- GammaArray[, z[rev_t+1], rev_t] * alpha[rev_t, ]
      trans_probs <- trans_probs / sum(trans_probs)

      # Sample state directly
      z[rev_t] <- rcat(1, trans_probs)
    }

    returnType(double(1))
    return(z)
  }
)

