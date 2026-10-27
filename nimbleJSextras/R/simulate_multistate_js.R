
#' @title Simulate Multistate Jolly-Seber Capture-Recapture Data
#' @param N Superpopulation size
#' @param xi Vector of entry probabilities, not necessarily normalized (length = `K`).
#' @param alpha A matrix of initial state assignment probabilites upon entry (dim = `J-2`x`K`). The total
#' number of states including pre-birth and death is `J`.
#' @param phi A `J-2`x`K-1` matrix of survival probabilities for each of the `J-2` alive states
#' @param psi A `J-2`x`J-2`x`K-1` array where `psi[,,t]` is the probability transition matrix between the `J-2`
#' alive states.
#' @param p A `J-2`x`O`x`K` array where `p[i,j,t]` is the probability of making observation `j` when the
#' individual is in true state `i` on occasion `t`. The number of observation states is `O`. The observation states
#' do not need to be directly related the the true states.
#' @import nimble
#' @export
simulate_multistate_js <- nimbleFunction(
  run = function(N = integer(0),
                 xi = double(1),     # K vector
                 alpha = double(2),  # (J-2) x K
                 phi = double(2),    # (J-2) x (K-1)
                 psi = double(3),    # (J-2) x (J-2) x (K-1)
                 p = double(3)       # (J-2) x O x K
  ) {

    returnType(sim_js_list())
    J <- dim(p)[1] + 2
    O <- dim(p)[2]
    K <- dim(p)[3]
    z <- matrix(1,N,K)
    x <- matrix(0,N,K)
    observed <- integer(N)

    # Conditional entry probabilities
    xi_tilde <- numeric(K)
    for(t in 1:(K-1)){xi_tilde[t] <- xi[t]/sum(xi[t:K])}
    if(xi[K]>0){
      xi_tilde[K] <- 1
    }else{
      xi_tilde[K] <- 0
    }

    # initial probabilities
    pi <- numeric(J)
    pi[1] <- 1-xi_tilde[1]
    pi[2:(J-1)] <- xi_tilde[1]*alpha[1:(J-2),1]
    pi[J] <- 0

    #TPM
    G <- array(0.0, dim=c(J,J,K-1))
     for(t in 1:(K-1)){
      G[1,1,t] <- 1-xi_tilde[t+1]
      for(l in 2:(J-1)){
        G[1,l,t] <- xi_tilde[t+1]*alpha[l-1, t+1]
        G[l, 2:(J-1), t] <- phi[l-1, t] * psi[l-1, 1:(J-2), t]
        G[l, J, t] <- 1-phi[l-1, t]
        }
      G[J,J,t] <- 1
    }

    #Sim states and ch
    for (i in 1:N) {
      z[i,1] <- rcat(1, pi)
      if((z[i,1]>1) & (z[i,1]<J)){
        x[i,1] <- rcat(1, p[z[i,1]-1, 1:O, 1])
      } else{
        x[i,1] <- 0
      }
      if(x[i,1]>0) observed[i] <- 1
      for(t in 2:K){
        z[i,t] <- rcat(1, G[z[i,t-1], 1:J, t-1])
        if((z[i,t]>1) & (z[i,t]<J)){
          x[i,t] <- rcat(1, p[z[i,t]-1, 1:O, t])
        } else{
          x[i,t] <- 0
        }
        if((observed[i]==1) | (x[i,t]>0)) observed[i] <- 1
      }
    }
    out <- sim_js_list$new()
    out$ch <- x
    out$true_states <- z
    out$observed <- observed

    return(out)
  }
)
