#' @export
#' @import nimble
load_Gamma <- nimbleFunction(
  run = function(xi_tilde = double(), phi = double()) {
    returnType(double(2))

    Gamma <- nimMatrix(0,3,3)
    Gamma[1,1] <- 1 - xi_tilde
    Gamma[1,2] <- xi_tilde
    Gamma[2,2] <- phi
    Gamma[2,3] <- 1- phi
    Gamma[3,3] <- 1
    return(Gamma)
  })
