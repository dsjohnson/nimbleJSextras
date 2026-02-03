#' @export
#' @import nimble
load_mix_Gamma <- nimbleFunction(
  run = function(xi_tilde = double(1), phi = double(1)) {
    returnType(double(2))

    nmix <- length(phi)
    idx <- c(1:(3*nmix))[rep(1:3, nmix)==1]-1
    Gamma <- nimMatrix(0,3*nmix,3*nmix)
    for(i in 1:nmix){
      Gamma[idx[i]+1,idx[i]+1] <- 1 - xi_tilde[i]
      Gamma[idx[i]+1,idx[i]+2] <- xi_tilde[i]
      Gamma[idx[i]+2,idx[i]+2] <- phi[i]
      Gamma[idx[i]+2,idx[i]+3] <- 1- phi[i]
      Gamma[idx[i]+3,idx[i]+3] <- 1
    }
    return(Gamma)
  })
