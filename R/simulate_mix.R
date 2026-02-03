#' @import nimble
#' @export
simList <- nimbleList(
  id_obs=integer(1), ch = double(2), states = double(2)
)


#' @export
#' @import nimble
simulate_mix_ch <- nimbleFunction(
  run = function(
    n = integer(default=1),
    alpha = double(1),
    xi = double(2),
    phi = double(2),
    p = double(2)
  ) {

    returnType(simList())
    res <- simList$new()

    K <- dim(xi)[1]
    nmix <- dim(xi)[2]
    ch <- nimMatrix(nrow=n, ncol=K)
    z <- nimMatrix(nrow=n, ncol=K)
    idx <- c(1:(3*nmix))[rep(1:3, nmix)==1]-1

    xi_tilde <- xi
    for(j in 1:nmix){
      for(t in 1:K){
        xi_tilde[t,j] <- xi_tilde[t,j]/sum(xi_tilde[t:K,j])
      }
    }

    pi <- rep(0,3*nmix)
    for(j in 1:nmix){
      pi[idx[j]+1] <- alpha[j]*(1-xi_tilde[1,j])
      pi[idx[j]+2] <- alpha[j]*xi_tilde[1,j]
    }

    Gamma <- nimArray(dim=c(3*nmix, 3*nmix, K-1))
    for(j in 1:nmix){
      for(t in 1:(K-1)){
        Gamma[idx[j]+1,idx[j]+1, t] <- 1 - xi_tilde[t+1,j]
        Gamma[idx[j]+1,idx[j]+2, t] <- xi_tilde[t+1,j]
        Gamma[idx[j]+2,idx[j]+2, t] <- phi[t,j]
        Gamma[idx[j]+2,idx[j]+3, t] <- 1- phi[t,j]
        Gamma[idx[j]+3,idx[j]+3, t] <- 1
      }
    }

    p_det <- nimMatrix(nrow=K, ncol=K)
    for(j in 1:nmix){p_det[,idx[j]+2] <- p[,j]}

    for(i in 1:n){
      z[i,1] <- rcat(1,pi)
      ch[i,1] <- rbinom(1, 1, p_det[1,z[i,1]])
      for(t in 2:K){
        z[i,t] <- rcat(1,Gamma[z[i,t-1],,t-1])
        ch[i,t] <- rbinom(1, 1, p_det[t,z[i,t]])
      }
    }
    det <- rep(0,n)
    for(i in 1:n){det[i] <- sum(ch[i,])>0}
    id_obs <- c(1:n)[det==1]
    res$id_obs <- id_obs
    res$ch <- ch[id_obs,]
    res$states <- z
    return(res)
  })
