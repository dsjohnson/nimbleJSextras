#' Sample latent states for undetected (or detected) individuals in Jolly-Seber models
#'
#' Forward-filtering backward-sampling functions to sample sequence of latent states
#' given capture histories for Jolly-Seber HMM models.
#'
#' @aliases sample_det sample_det_ms sample_det_binom sample_det_pois sample_det_mark
#' @name sample_det
#'
#' @param x Vector of observed counts or capture status for an individual across time (length `K`).
#' @param init Vector of initial state probabilities at time `t = 1` of length `J+2`
#'   (where `J` is the number of intermediate active states, plus 1 initial unentered and 1 final exited state).
#' @param probObs Time-dependent 3D array of observation probabilities. First
#'   two dimensions of `probObs` are of size (number of possible system states)
#'   `x` (number of possible observation classes) and the third dimension is the
#'   number of occasions (`K`).
#' @param probTrans Time-dependent 3D array of system state transition probabilities
#'   with dimensions `(J+2) x(J+2) x K` (where `J+2` is the total number of system states and `K` is number of occasions).
#'   `probTrans[r, s, t]` is the probability of transitioning from state `r` at time `t` to state `s` at time `t+1`.
#' @param prob 2D array (`K x J`) containing detection probabilities for each
#'   intermediate/active state `J` at time `t`.
#' @param size Vector of sizes (number of trials) for the binomial observation at each time (length `K`).
#' @param rate 2D array (`K x J`) containing Poisson detection rates for each
#'   intermediate/active state `J` at time `t`.
#' @param mark Vector of continuous observations/marks across time (length `K`).
#' @param knots Vector of knot locations/means corresponding to each intermediate state `J` (length `J`).
#' @param sigma_knot Standard deviation around knot locations.
#'
#' @return A numeric vector of length `K` containing the sampled latent state sequence.
#'
#' -----------------------------------------------------------------------------
#' Multistate JS
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname sample_det
sample_det_ms <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 probObs = double(3),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Dimensions
    n_states <- length(init)
    K <- length(x)  # Ensure len matches obs length
    # Forward probabilities
    alpha <- matrix(0, nrow=K, ncol=n_states)
    sampled_states <- numeric(K)  # To store sampled states
    # ---- Forward Pass ----
    for (t in 1:K) {
      if(t == 1){
        alpha[t, ] <- init * probObs[, x[t], t]
      } else{
        alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * probObs[, x[t], t]
      }
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }
    # ---- Backward Sampling ----
    sampled_states[K] <- rcat(1, alpha[K, ])  # Sample final state
    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional backward probabilities
      back_probs <- probTrans[, sampled_states[rev_t+1], rev_t] * alpha[rev_t, ]
      back_probs <- back_probs / sum(back_probs)
      # Sample state directly
      sampled_states[rev_t] <- rcat(1, back_probs)
    }
    returnType(double(1))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' -----------------------------------------------------------------------------
#' Binomial obs (robust design) JS
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @importFrom stats dnorm
#' @rdname sample_det
sample_det_binom <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 prob = double(2),     # Emission probabilities
                 size = double(1),
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Dimensions
    n_states <- length(init)
    K <- length(x)  # Ensure len matches obs length
    J <- dim(prob)[2]
    # Forward probabilities
    alpha <- matrix(0, nrow=K, ncol=n_states)
    sampled_states <- numeric(K)  # To store sampled states
    Pdiag <- nimNumeric(value = 0.0, length = n_states)
    # ---- Forward Pass ----
    for (t in 1:K) {
      xt <- x[t]
      Rt <- size[t]
      # if (xt == 0) {
      #   Pdiag[1] <- 1.0
      #   Pdiag[n_states] <- 1.0
      # } else {
      #   Pdiag[1] <- 0.0
      #   Pdiag[n_states] <- 0.0
      # }
      # ln_binom_coef <- 0.0
      # if (xt > 0 & xt < Rt) {
      #   ln_binom_coef <- lgamma(Rt + 1.0) - lgamma(xt + 1.0) - lgamma(Rt - xt + 1.0)
      # } else{
      #   ln_binom_coef <- 0.0
      # }
      for (k in 1:n_states) {
        # p_jt <- prob[t, j]
        # ln_dens <- ln_binom_coef + xt * log(p_jt) + (Rt - xt) * log(1.0 - p_jt)
        if(k==1 | k==n_states){
          Pdiag[k] <-  dbinom(xt, Rt, 0.0)
        } else {
          Pdiag[k] <-  dbinom(xt, Rt, prob[t,k-1])
        }
      }
      # for(k in 1:J) Pdiag[k+1] <- dbinom(x[t], size[t], prob[t,k], log=0)
      if(t ==1){
        alpha[t, ] <- init * Pdiag
      } else {
        alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      }
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }
    # ---- Backward Sampling ----
    sampled_states[K] <- rcat(1, alpha[K, ])  # Sample final state
    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      back_probs <- probTrans[, sampled_states[rev_t+1], rev_t] * alpha[rev_t, ]
      back_probs <- back_probs / sum(back_probs)
      # Sample state directly
      sampled_states[rev_t] <- rcat(1, back_probs)
    }
    returnType(double(1))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' -----------------------------------------------------------------------------
#' Poisson JS (continuous-time robust design)
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname sample_det
sample_det_pois <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 init = double(1),        # Initial state probabilities
                 rate = double(2),     # Emission probabilities
                 probTrans = double(3)   # Transition probabilities
  ) {
    # Dimensions
    n_states <- length(init)
    K <- length(x)  # Ensure len matches obs length
    J <- dim(rate)[2]
    # Forward (alpha) matrix
    alpha <- matrix(0, nrow=K, ncol=n_states)
    sampled_states <- numeric(K)  # To store sampled states
    lambda <- nimNumeric(value = 0.0, length = n_states)
    # ---- Forward Pass ----
    for (t in 1:K) {
      lambda[2:(J+1)] <- rate[t,]
      Pdiag <- dpois(x[t], lambda, log=0)
      if(t==1){
        alpha[t, ] <- init * Pdiag
      } else{
        alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      }
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }
    # ---- Backward Sampling ----
    sampled_states[K] <- rcat(1, alpha[K, ])  # Sample final state
    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      back_probs <- probTrans[, sampled_states[rev_t+1], rev_t] * alpha[rev_t, ]
      back_probs <- back_probs / sum(back_probs)

      # Sample state directly
      sampled_states[rev_t] <- rcat(1, back_probs)
    }
    returnType(double(1))  # Return sampled state sequence
    return(sampled_states)
  }
)

#' -----------------------------------------------------------------------------
#' Bernoulli/Gaussian marked obs JS
#' -----------------------------------------------------------------------------
#' @export
#' @import nimble
#' @rdname sample_det
sample_det_mark <- nimble::nimbleFunction(
  run = function(x = double(1),         # Observation sequence
                 mark = double(1),
                 init = double(1),        # Initial state probabilities
                 prob = double(2),     # Emission probabilities
                 probTrans = double(3),   # Transition probabilities
                 knots = double(1),
                 sigma_knot = double()
  ) {
    returnType(double(1))  # Return sampled state sequence
    # Dimensions
    J <- length(knots)
    n_states <- J+2
    K <- length(x)  # Ensure len matches obs length
    # Forward probabilities
    alpha <- matrix(0, nrow=K, ncol=n_states)
    sampled_states <- numeric(K)  # To store sampled states
    dm <- nimNumeric(value = 1.0, length = n_states)
    pt <- nimNumeric(value = 0.0, length = n_states)
    # ---- Forward Pass ----
    ### Observation distribution -----------
    if(!is.na(mark[1])){
      dm[2:(J+1)] <- dnorm(mark[1], mean=knots, sd=sigma_knot, log=FALSE)
    } else {
      dm[2:(J+1)] <- 1.0
    }
    pt[2:(J+1)] <- prob[1,]
    Pdiag <- (1-x[1])*(1-pt) + x[1]*pt*dm
    ### ------------------------------------
    alpha[1, ] <- init * Pdiag
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
    for (t in 1:K) {
      ### Observation distribution -----------
      if(!is.na(mark[t])){
        dm[2:(J+1)] <- dnorm(mark[t], mean=knots, sd=sigma_knot, log=FALSE)
      } else {
        dm[2:(J+1)] <- 1.0
      }
      pt[2:(J+1)] <- prob[t,]
      Pdiag <- (1-x[t])*(1-pt) + x[t]*pt*dm
      ### ------------------------------------
      if(t==1){
        alpha[t, ] <- init * Pdiag
      } else{
        alpha[t, ] <- (alpha[t-1, ] %*% probTrans[, ,t-1]) * Pdiag
      }
      alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
    }
    # ---- Backward Sampling ----
    sampled_states[K] <- rcat(1, alpha[K, ])  # Sample final state
    for (t in 1:(K-1)) {
      rev_t <- K-t # rev_t goes from K-1 to 1
      # Compute conditional probabilities
      back_probs <- probTrans[, sampled_states[rev_t+1], rev_t] * alpha[rev_t, ]
      back_probs <- back_probs / sum(back_probs)
      # Sample state directly
      sampled_states[rev_t] <- rcat(1, back_probs)
    }
    return(sampled_states)
  }
)
