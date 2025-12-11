iar_Q <- function(n,p=1)
{
  if(n<2*p) stop("n must be >= 2*p\n")
  tmp1 <- out <- matrix(0,n,n)
  tmp2 <- ((-1)^c(0:p))*choose(p,c(0:p))
  for(i in (p+1):n)
  {
    tmp1[,i] <- c(rep(0,i-p-1), tmp2, rep(0, n-i))
  }
  for(i in n:(p+1))
  {
    tmp4 <- tmp1[,c(1:n)[tmp1[i,]!=0]]
    tmp5 <- tmp1[i,tmp1[i,]!=0]
    tmp6 <- t(t(tmp4) * tmp5)
    tmp6[i,] <- 0
    out[i,] <- -rowSums(tmp6)
  }
  out[1:p,] <- out[n:(n-p+1),n:1]
  return(diag(apply(out,1,sum)) - out)
}


Q_to_car <- function(Q, tol = 1e-10) {
  # Ensure Q is a matrix
  Q <- as.matrix(Q)

  # Validation
  if (nrow(Q) != ncol(Q)) stop("Q must be a square matrix.")
  if (!isSymmetric(Q, tol = tol)) warning("Q is not symmetric. dcar_normal requires a symmetric structure.")

  N <- nrow(Q)
  adj <- numeric()
  weights <- numeric()
  num <- integer(N)

  for (i in 1:N) {
    # 1. Identify Neighbors
    # Find column indices where Q[i,j] is non-zero, excluding the diagonal (i)
    # Using 'tol' to handle floating point zeros
    nbs <- which(abs(Q[i, ]) > tol)
    nbs <- nbs[nbs != i] # Remove self-loop

    # 2. Store Number of Neighbors
    num[i] <- length(nbs)

    # 3. Store Adjacency and Weights
    if (length(nbs) > 0) {
      adj <- c(adj, nbs)

      # IMPORTANT: nimble weights are -Q[i, j]
      # The diagonal Q[i,i] is ignored here because nimble calculates it internally
      # as the sum of the weights.
      w_vals <- -Q[i, nbs]
      weights <- c(weights, w_vals)
    }
  }

  return(list(
    adj = adj,
    weights = weights,
    num = num
  ))
}

#' @title Create liest arguments for a random walk of order `p`.to be used
#' in the `dcar_normal` distribution
#'
#' @description This function calculates the entries for a precision matrix for a random walk of order `p`.
#'
#' @param n The length of the RW(p) process. The length must be greater than or equal to \code{2*p}.
#' @param p The order of the random walk process.
#'
#' @return A list with elements: `adj`, `weights`, and `num` for use with the
#'  `dcar_normal` distribution in `nimble` models
#'
#' @references H. Rue and L. Held (2005) Gaussian Markov Random Fields. Chapman & Hall/CRC. 263 pp.
#' @export

make_nimble_icar <- function(n,p=1){
  Q <- iar_Q(n,p)
  out <- Q_to_car(Q)
  out$c <- p
  return(out)
}
