#' @export
#' @import nimble
load_mix_Gamma <- nimbleFunction(
  run = function(xi_tilde = double(1), phi = double(1)) {
    returnType(double(2))

    nmix <- length(phi)
    idx <- c(1:nmix^2)[rep(1:3, 3)==1]-1
    Gamma <- nimMatrix(0,nmix^2,nmix^2)
    for(i in 1:nmix){
      Gamma[idx[i]+1,idx[i]+1] <- 1 - x_tilde[i]
      Gamma[idx[i]+1,idx[i]+2] <- x_tilde[i]
      Gamma[idx[i]+2,idx[i]+2] <- phi[i]
      Gamma[idx[i]+2,idx[i]+3] <- 1- phi[i]
      Gamma[idx[i]+3,idx[i]+3] <- 1
    }
    return(Gamma)
  })
