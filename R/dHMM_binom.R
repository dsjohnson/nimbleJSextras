#' @export
#' @rdname dDHMM_bern
dDHMMo_bern <- nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 len = integer(0, default=0),## length of x (needed as a separate param for rDHMMo_binom)
                 log = integer(0, default = 0)) {
    pi <- init # State probabilities at time t=1
    logL <- 0
    len <- length(x)
    for (t in 1:len) {
      ### Observation distribution -----------
      Pdiag <- dbinom(x[t], 1, prob[t,])
      ### ------------------------------------
      Zpi <- Pdiag * pi
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t != len) pi <- ((Zpi %*% probTrans[,,t]) / sumZpi)[1, ]
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
#' @rdname dDHMM_bern
rDHMMo_bern <- nimbleFunction(
  run = function(n = integer(),
                 init = double(1),
                 prob = double(2),
                 probTrans = double(3),
                 len = integer(0, default=0),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    x <- numeric(len)
    state <- rcat(1,init)
    for(t in 1:len){
      x[t] <- rbinom(1,1,prob[t,state])
      if(t!=len) state <- rcat(1,probTrans[state,,t])
    }
    returnType(double(1))
    return(x)
  }
)


#' @export
#' @rdname dDHMM_binom
dDHMMo_binom <- nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 len = integer(0, default=0),## length of x (needed as a separate param for rDHMMo_binom)
                 log = integer(0, default = 0)) {
    pi <- init # State probabilities at time t=1
    logL <- 0
    len <- length(x)
    for (t in 1:len) {
      ### Observation distribution -----------
      Pdiag <- dbinom(x[t], size[t], prob[t,])
      ### ------------------------------------
      Zpi <- Pdiag * pi
      sumZpi <- sum(Zpi)
      logL <- logL + log(sumZpi)
      if (t != len) pi <- ((Zpi %*% probTrans[,,t]) / sumZpi)[1, ]
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
#' @rdname dDHMM_binom
rDHMMo_binom <- nimbleFunction(
  run = function(n = integer(),
                 init = double(1),
                 prob = double(2),
                 size = double(1),
                 probTrans = double(3),
                 len = integer(0, default=0),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    x <- numeric(len)
    state <- rcat(1,init)
    for(t in 1:len){
      x[t] <- rbinom(1,size[t],prob[t,state])
      if(t!=len) state <- rcat(1,probTrans[state,,t])
    }
    returnType(double(1))
    return(x)
  }
)
