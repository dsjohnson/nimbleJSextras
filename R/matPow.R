#' @export
#' @import nimble
matPow <- nimbleFunction(
  run = function(A = double(2), p = integer(0)) {
    returnType(double(2))

    # Square matrix check (handled by compiler if dimensions are static)
    n <- dim(A)[1]
    res <- diag(n)
    base <- A

    # Binary exponentiation loop
    # Exponent is reduced by half each step: O(log p) multiplications
    while(p > 0) {
      if(p %% 2 == 1) {
        res <- res %*% base
      }
      base <- base %*% base
      p <- floor(p / 2)
    }

    return(res)
  }
)
