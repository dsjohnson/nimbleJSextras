#' @export
#' @import nimble
load_mix_P <- nimbleFunction(
  run = function(p = double(1)) {
    returnType(double(2))

    nmix <- length(p)
    idx <- c(1:(3*nmix))[rep(1:3, nmix)==1]-1
    P <- nimMatrix(0,3*nmix,2)
    P[,1] <- 1
    for(i in 1:nmix){
      P[idx[i]+2,1] <- 1 - p[i]
      P[idx[i]+2,2] <- p[i]
    }
    return(P)
  })
