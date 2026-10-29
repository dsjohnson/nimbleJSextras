
predict_abundance <- nimbleFunction(
  run = function(xi=double(2), phi = double(1), tau=double(2), p=double(2),
                 eta=double(2), Nu = double(1), x = double(2)) {
    iter <- dim(xi)[1]
    K <- dim(xi)[2]
    n <- dim(x)[1]
    m <- dim(tau)[2]
    n_states <- 6+2*m

    xi_tilde <- numeric(K)
    Gamma <- nimArray(dim=c(n_states,n_states,K-1))
    pi <- numeric(n_states)
    Pmats <- nimArray(dim=c(n_states,3,K))
    Nu_t <- nimMatrix(nrow=n_states, ncol=K)
    z <- nimMatrix(nrow=n, ncol = K)
    n_t <- numeric(K)
    N_t <- nimMatrix(nrow=iter, ncol=K)

    for(j in 1:iter){

      ### Load HMM matrices
      for(t in 1:(K-1)){
        xi_tilde[t] <- xi[j,t]/sum(xi[j,t:K])
      }
      xi_tilde[K] <- 1


      for(t in 1:(K-1)){
        Gamma[,,t] <- make_honu_Gamma(xi_tilde[t+1], phi[j], p[j,t+1], eta[j,t+1], tau[j,])
        Pmats[,,t] <- make_honu_P(m)
      }
      Pmats[,,K] <- make_honu_P(m)
      pi <- make_honu_pi(xi_tilde[1], eta[j,1], m)

      ### Predict abundance for uncaptured id
      Nu_t <- rJSUndet_ms(1, Nu[j], pi, Pmats, Gamma)

      ### Predict state for captured id
      for(i in 1:n){
        z[i, ] <- rJSCond_ms(1, x[i, ], pi, Pmats, Gamma)
      }
      ### Determine abundance 
      for(t in 1:K){
        n_t[t] <- n - (sum(z[,t]==1) + sum(z[,t]==n_states))
        N_t[j, t] <-  n_t[t] + (Nu[j] - (Nu_t[1, t] + Nu_t[n_states, t]))
      }

    } # MCMC iteration loop end (j index)

    returnType(double(2))
    return(N_t)
  }) # function end
