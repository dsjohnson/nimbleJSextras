# #'Convenience function for loading parameter arrays for JS models
# #'@aliases load_Gamma load_mix_Gamma load_mix_P
# #'@name load_foo
#'
# #'@export
# #'@import nimble
# load_Gamma <- nimbleFunction(
#   run = function(xi_tilde = double(), phi = double()) {
#     returnType(double(2))
#
#     Gamma <- nimMatrix(0,3,3)
#     Gamma[1,1] <- 1 - xi_tilde
#     Gamma[1,2] <- xi_tilde
#     Gamma[2,2] <- phi
#     Gamma[2,3] <- 1- phi
#     Gamma[3,3] <- 1
#     return(Gamma)
#   })


# #'@export
# #'@import nimble
# #'@rdname load_foo
# load_mix_Gamma <- nimbleFunction(
#   run = function(xi_tilde = double(1), phi = double(1)) {
#     returnType(double(2))
#
#     nmix <- length(phi)
#     idx <- c(1:(3*nmix))[rep(1:3, nmix)==1]-1
#     Gamma <- nimMatrix(0,3*nmix,3*nmix)
#     for(i in 1:nmix){
#       Gamma[idx[i]+1,idx[i]+1] <- 1 - xi_tilde[i]
#       Gamma[idx[i]+1,idx[i]+2] <- xi_tilde[i]
#       Gamma[idx[i]+2,idx[i]+2] <- phi[i]
#       Gamma[idx[i]+2,idx[i]+3] <- 1- phi[i]
#       Gamma[idx[i]+3,idx[i]+3] <- 1
#     }
#     return(Gamma)
#   })


# #'@export
# #'@import nimble
# #'@rdname load_foo
# load_mix_P <- nimbleFunction(
#   run = function(p = double(1)) {
#     returnType(double(2))
#
#     nmix <- length(p)
#     idx <- c(1:(3*nmix))[rep(1:3, nmix)==1]-1
#     P <- nimMatrix(0,3*nmix,2)
#     P[,1] <- 1
#     for(i in 1:nmix){
#       P[idx[i]+2,1] <- 1 - p[i]
#       P[idx[i]+2,2] <- p[i]
#     }
#     return(P)
#   })
