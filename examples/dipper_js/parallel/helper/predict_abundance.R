

predict_abundance <- nimbleFunction(
  run = function(xi=double(2), phi = double(2), p=double(2),
                 Nu = double(1), x = double(2)) {
    iter <- dim(xi)[1]
    K <- dim(xi)[2]
    n <- dim(x)[1]

    xi_tilde <- numeric(K)
    Gamma <- nimArray(dim=c(3,3,K-1))
    pi <- numeric(3)
    Pmats <- nimArray(dim=c(3,3,K))
    Nu_t <- nimMatrix(nrow=3, ncol=K)
    Nu_t_active <- numeric(K)
    z <- nimMatrix(nrow=n, ncol = K)
    N_t <- nimMatrix(nrow=iter, ncol=K)

    for(j in 1:iter){

      ### Load HMM matrices
      for(t in 1:(K-1)){
        xi_tilde[t] <- xi[j,t]/sum(xi[j,t:K])
      }
      xi_tilde[K] <- 1

      pi[1] <- 1-xi_tilde[1]
      pi[2] <- xi_tilde[1]
      pi[3] <- 0

      for(t in 1:K){Pmats[1:3,1:2,t] <- load_probObs(p[j,t])}

      for(t in 1:(K-1)){
        Gamma[1:3, 1:3, t] <- load_Gamma(xi_tilde[t+1], phi[j,t])
      }

      ### Predict abundance for uncaptured id
      Nu_t <- rJSUndet_ms(1, Nu[j], pi, Pmats, Gamma)
      Nu_t_active <- Nu_t[2,]

      ### Predict state for captured id
      for(i in 1:n){
        z[i, ] <- rJSCond_ms(1, x[i, ], pi, Pmats, Gamma)
      }
      ### Determine abundance for captured animals
      for(t in 1:K){
        N_t[j, t] <-  sum(z[,t]==2) + Nu_t_active[t]
      }

    } # MCMC iteration loop end (j index)
    returnType(double(2))
    return(N_t)
  }) # function end
