

predict_abundance <- nimbleFunction(
  run = function(xi=double(2), phi = double(2), p=double(2),
                 nu = double(1), x = double(2)) {
    iter <- dim(xi)[1]
    K <- dim(xi)[2]
    n <- dim(x)[1]

    xi_tilde <- numeric(K)
    Gamma <- nimArray(dim=c(3,3,K-1))
    pi <- numeric(3)
    Pmats <- nimArray(dim=c(3,3,K))
    nu_t <- nimMatrix(nrow=3, ncol=K)
    n_avail_u <- numeric(K)
    z <- nimMatrix(nrow=n, ncol = K)
    avail_obs <- nimMatrix(nrow=n, ncol = K)
    N_avail <- nimMatrix(nrow=iter, ncol=K)

    for(j in 1:iter){

      ### Load HMM matrices
      for(t in 1:(K-1)){
        xi_tilde[t] <- xi[j,t]/sum(xi[j,t:K])
      }
      xi_tilde[K] <- 1

      pi[1] <- 1-xi_tilde[1]
      pi[2] <- xi_tilde[1]
      pi[3] <- 0

      for(t in 1:K){
        Pmats[1,1,t] <- 1
        Pmats[1,2,t] <- 0
        Pmats[2,1,t] <- 1-p[j,t]
        Pmats[2,2,t] <- p[j,t]
        Pmats[3,1,t] <- 1
        Pmats[3,2,t] <- 0
      }

      for(t in 1:(K-1)){
        Gamma[1,1,t] <- 1-xi_tilde[t+1]
        Gamma[1,2,t] <- xi_tilde[t+1]
        Gamma[1,3,t] <- 0
        Gamma[2,1,t] <- 0
        Gamma[2,2,t] <- phi[j,t]
        Gamma[2,3,t] <- 1 - phi[j,t]
        Gamma[3,1,t] <- 0
        Gamma[3,2,t] <- 0
        Gamma[3,3,t] <- 1
      }

      ### Predict abundance for uncaptured id
      nu_t <- sample_undet_ms(nu[j], pi, Pmats, Gamma)
      n_avail_u <- nu_t[2,]

      ### Predict state for captured id
      for(i in 1:n){
        z[i, ] <- sample_det_ms(x[i, ], pi, Pmats, Gamma)

      }
      ### Determine abundance for captured animals
      for(t in 1:K){
        n_avail_obs <- 0
        for(i in 1:n){
          # Check if individual is 'available' (state = 2)
          if(z[i,t]==2) {n_avail_obs <- n_avail_obs + 1}
        }
        # Assign to your result matrices
        N_avail[j, t] <-  n_avail_obs + n_avail_u[t]
      }

    } # MCMC iteration loop end (j index)
    returnType(double(2))
    return(N_avail)
  }) # function end
