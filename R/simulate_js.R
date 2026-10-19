#' @title Simulate Traditional Jolly-Seber Capture-Recapture Data
#' @param N Superpopulation size
#' @param xi Vector of entry probabilities, not necessarily normalized (length = `K`).
#' @param phi A `K-1` vector of survival probabilities.
#' @param p A `K` vector of capture probabilities for when an individual is available for capture.
#' @import nimble
#' @export
simulate_js <- nimbleFunction(
  run = function(N = integer(0),
                 xi = double(1),     # K vector
                 phi = double(1),    # K-1
                 p = double(1)       # K
  ) {

    returnType(sim_js_list())
    K <- length(p)
    z <- matrix(1,N,K)
    x <- matrix(0,N,K)
    observed <- integer(N)

    # conditional entry probabilities
    xi_tilde <- numeric(K)
    for(t in 1:(K-1)){xi_tilde[t] <- xi[t]/sum(xi[t:K])}
    if(xi[K]>0){
      xi_tilde[K] <- 1
    }else{
      xi_tilde[K] <- 0
    }

    # initial probabilities
    pi <- numeric(3)
    pi[1] <- 1-xi_tilde[1]
    pi[2] <- xi_tilde[1]
    pi[3] <- 0

    #TPM
    # G <- array(0.0, dim=c(3,3,K-1))
    # for(t in 1:(K-1)){
    #   G[1,1,t] <- 1-xi_tilde[t+1]
    #   G[1,2,t] <- xi_tilde[t+1]
    #   G[2, 2, t] <- phi[t]
    #   G[2, 3, t] <- 1-phi[t]
    #   G[3,3,t] <- 1
    # }

    #Sim states and ch
    for (i in 1:N) {
      # sample state for idiv i
      z[i,1] <- rcat(1, pi)
      if(z[i,1]==2){
        x[i,1] <- rbinom(1, 1, p[1])
      } else{
        x[i,1] <- 0
      }
      if(x[i,1]==1) observed[i] <- 1
      for(t in 2:K){
        # sample state
        if(z[i,t-1]==1){
          z[i,t] <- rbinom(1, 1,  xi_tilde[t]) + 1
        } else if(z[i,t-1]==2){
          z[i,t] <- z[i,t-1] + (1-rbinom(1,1,phi[i,t-1]))
        } else{
          z[i,t] <- z[i,t-1]
        }
        # sample obs for indiv i
        if(z[i,t]==2){
          x[i,t] <- rbinom(1, 1, p[t])
        } else{
          x[i,t] <- 0
        }
        if((observed[i]==1) | (x[i,t]==1)) observed[i] <- 1
      }
    }
    out <- sim_js_list$new()
    out$ch <- x
    out$true_states <- z
    out$observed <- observed

    return(out)
  }
)
