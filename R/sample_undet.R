
#' @import nimble
#' @export
sample_state_undet <- nimble::nimbleFunction(
  run = function(
    nu = integer(),           # Total number of undetected indiv.
    len = integer(),
    init = double(1),        # Initial state probabilities
    probTrans = double(2)    # Transition probabilities
  ) {
    J <- dim(probTrans)[1]
    K <- len
    nu_t <- matrix(0,nrow=J, ncol=K)
    Tu <- array(0,dim = c(J,J,K-1))
    nu_t[1:J,1] <- rmulti(1, nu, init[1:J])
    for(t in 2:K){
      for(l in 1:J){
        Tu[l,1:J,t-1] <- rmulti(1,  nu_t[l,t-1], probTrans[l,1:J])
        nu_t[l,t] <- sum(Tu[1:J,l,t-1])
      }
    }
    returnType(double(2))  # Return sampled state sequence
    return(nu_t)
  }
)



#' @import nimble
#' @export
sample_state_undet_D <- nimble::nimbleFunction(
  run = function(
    nu = integer(),           # Total number of undetected indiv.
    init = double(1),        # Initial state probabilities
    probTrans = double(3)    # Transition probabilities
  ) {
    J <- dim(probTrans)[1]
    K <- dim(probTrans)[3] + 1
    nu_t <- matrix(0,nrow=J, ncol=K)
    Tu <- array(0,dim = c(J,J,K-1))
    nu_t[1:J,1] <- rmulti(1, nu, init[1:J])
    for(t in 2:K){
      for(l in 1:J){
        Tu[l,1:J,t-1] <- rmulti(1,  nu_t[l,t-1], probTrans[l,1:J,t-1])
        nu_t[l,t] <- sum(Tu[1:J,l,t-1])
      }
    }
    returnType(double(2))  # Return sampled state sequence
    return(nu_t)
  }
)
