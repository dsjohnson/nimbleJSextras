predict_abundance <- nimbleFunction(
  run = function(xi=double(2), phi = double(2), theta=double(1), p=double(2),
                 nd=integer(), nu = integer(1), x = integer(2)) {
    iter <- dim(xi)[1]
    K <- dim(xi)[2]
    n <- dim(x)[1]
    n_states <- 3+ nd

    x_tilde <- numeric(K)
    Gamma <- nimArray(dim=c(n_states,n_states,K-1))
    pi <- numeric(n_states)
    Pmats <- nimArray(dim=c(n_states,2,K))


    for(j in 1:iter){

      ### Load HMM matrices
      for(t in 1:(K-1)){
        xi_tilde[t] <- xi[j,t]/sum(xi[j,t:K])

      }
      xi_tilde[K] <- 1

      ### Predict abundance

      nu_t[1:n_states, 1:K] <- sample_undet_ms(n=nu[j], init=pi,
                                               probObs=Pmats, probTrans=Gamma)
      for(i in 1:n){
        det_state[i,1:K] <- sample_det_ms(x=x, init=pi,
                                          probObs=Pmats, probTrans=Gamma)
      }
      nest_obs <- det_state==2
      avail_obs[1:nobs,1:K] <- 1-((det_state[1:nobs,1:K]==1) + (det_state[1:nobs,1:K]==nd))
      for(t in 1:K){
        N_nest[t] <- sum(nest_obs[1:nobs, t]) + nu_t[2,t]
        N[t] <- nu - (nu_t[1,t]+nu_t[nd,t]) + sum(avail_obs[1:nobs, t])
      }
    } # i loop end
  } # j loop end
