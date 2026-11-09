
matern_cor <- function(d, nu, rho) {
  if (d == 0) return(1)
  scaled_d <- sqrt(2 * nu) * d / rho
  const <- (2^(1 - nu)) / gamma(nu)
  const * (scaled_d^nu) * besselK(scaled_d, nu = nu)
}

#' Calibrate rho for Any nu and Knot Spacing
#'
#' @param delta Numeric. Spacing between adjacent knots.
#' @param nu Numeric. Smoothness parameter (nu > 0).
#' @param target_midpoint_cor Numeric. Target correlation at d = 0.5 * delta.
#'   Default is 0.75 for smooth overlap without severe collinearity.
#'
#' @return The calibrated value of rho.
#' @importFrom stats uniroot
#' @export
calibrate_matern_rho <- function(delta, nu = 1.5, target_midpoint_cor = 0.75) {
  midpoint_dist <- 0.5 * delta
  
  # Objective function: difference from target midpoint correlation
  objective <- function(rho) {
    matern_cor(d = midpoint_dist, nu = nu, rho = rho) - target_midpoint_cor
  }
  
  # Search within a safe multiplier range: [0.2 * delta, 10 * delta]
  fit <- uniroot(
    f = objective, 
    interval = c(0.2 * delta, 10 * delta), 
    tol = 1e-8
  )
  
  return(fit$root)
}


find_matern_padding <- function(nu, rho, threshold = 0.02) {

  # Find distance d where correlation drops to threshold
  uniroot(function(d) matern_cor(d, nu, rho) - threshold, 
          interval = c(0.1 * rho, 15 * rho))$root
}

#' @title Create Padded Knot Locations for Radial Basis Function Expansions
#'
#' @description Generates an equally spaced sequence of knot centers that extends 
#' beyond the boundary of the observed data domain \code{[min(x), max(x)]}. 
#' Padded "buffer" knots ensure overlapping basis coverage at the domain edges, 
#' preventing artificial collapse of prior variance in low-rank spline and 
#' Gaussian process models. Buffer width is dynamically determined based on the 
#' effective range of the Matérn kernel.
#'
#' @param x Numeric vector of observed data points or evaluation locations.
#' @param n_interior Integer. Number of knot locations desired within the interior 
#'   domain \code{[min(x), max(x)]}. Default is 15.
#' @param nu Numeric. Smoothness parameter of the Matérn basis kernel (\eqn{\nu > 0}). 
#'   Default is 1.5.
#' @param rho Numeric or \code{NULL}. Lengthscale/range parameter (\eqn{\rho > 0}). 
#'   If \code{NULL} (default), it is automatically calibrated using 
#'   \code{\link{calibrate_matern_rho}} to achieve an approximate partition of unity.
#' @param target_midpoint_cor Numeric. Target correlation at the midpoint between adjacent 
#'   knots (\eqn{0.5 \times \Delta}) used when calibrating \code{rho}. Default is 0.75. 
#'   Ignored if \code{rho} is supplied.
#' @param threshold Numeric. Correlation threshold defining the boundary cutoff distance.
#'  Default is 0.02 (0.02 correlation).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{knots}: Numeric vector of all knot locations, including interior and boundary buffer knots.
#'   \item \code{rho}: The calibrated or provided lengthscale/range parameter.
#'   \item \code{delta}: The spacing between adjacent knots.
#'   \item \code{n_pad}: Number of buffer knots added to each boundary.
#'   \item \code{pad_dist}: Distance extended beyond each boundary.
#' }
#'
#' @author Devin Johnson
#' @export
make_padded_knots <- function(x, 
                              n_interior = 15, 
                              nu = 1.5, 
                              rho = NULL, 
                              target_midpoint_cor = 0.75, 
                              threshold = 0.02) {
  x_min <- min(x)
  x_max <- max(x)
  
  # Step 1: Nominal knot spacing for the desired interior resolution
  delta <- (x_max - x_min) / (n_interior - 1)
  
  # Step 2: Calibrate rho if not explicitly provided
  if (is.null(rho)) {
    rho <- calibrate_matern_rho(delta = delta, nu = nu, target_midpoint_cor = target_midpoint_cor)
  }
  
  # Step 3: Compute exact padding distance based on tail decay threshold
  pad_dist <- find_matern_padding(nu = nu, rho = rho, threshold = threshold)
  n_pad <- ceiling(pad_dist / delta)
  
  # Step 4: Expanded knot sequence
  knots <- seq(x_min - n_pad * delta, x_max + n_pad * delta, by = delta)
  
  return(list(
    knots = knots, 
    rho = rho, 
    delta = delta, 
    n_pad = n_pad, 
    pad_dist = pad_dist
  ))
}



#' Matérn Radial Basis Function Matrix
#'
#' @description Generates a Matérn basis matrix \eqn{\Phi} for an input vector \code{x} 
#' evaluated across a set of knot locations \code{knots}.
#'
#' @param x Numeric vector of evaluation points (length \eqn{N}).
#' @param knots Numeric vector of basis knot centers (length \eqn{K}). 
#'   Defaults to an equally spaced sequence spanning \code{range(x)}.
#' @param nu Numeric. Smoothness parameter \eqn{\nu > 0}. Common choices are 
#'   \code{0.5} (Exponential / \eqn{C^0}), \code{1.5} (\eqn{C^1}), and \code{2.5} (\eqn{C^2}). 
#'   Default is \code{1.5}.
#' @param rho Numeric. Spatial range/lengthscale parameter \eqn{\rho > 0}. 
#'   Default is to optimize for the minimim size to provide a smooth function based
#'   on the value of `nu` selected.
#'
#' @return An \eqn{N \times K} numeric matrix where entry \eqn{(i, j)} is 
#'   \eqn{\text{Matérn}(|x_i - \text{knots}_j|; \nu, \rho)}.
#'
#' @export
matern_basis <- function(x, 
                         knots = seq(min(x), max(x), length.out = 10), 
                         nu = 1.5, 
                         rho = NULL) {
  
  if (nu <= 0) stop("nu must be > 0")
  if(!is.null(rho)){
    if (rho <= 0) stop("rho must be > 0")
  } else{
    delta <- min(diff(knots))
    rho <- calibrate_matern_rho(delta, nu = nu, target_midpoint_cor = 0.75)
  }
  
  # Absolute pairwise distance matrix: |x_i - knots_j| (N x K)
  d <- abs(outer(x, knots, "-"))
  
  # Closed-form implementations for common smoothness orders
  if (isTRUE(all.equal(nu, 0.5))) {
    # Exponential / Ornstein-Uhlenbeck (C^0)
    Phi <- exp(-d / rho)
    
  } else if (isTRUE(all.equal(nu, 1.5))) {
    # Matérn 3/2 (C^1)
    scaled_d <- sqrt(3) * d / rho
    Phi <- (1 + scaled_d) * exp(-scaled_d)
    
  } else if (isTRUE(all.equal(nu, 2.5))) {
    # Matérn 5/2 (C^2)
    scaled_d <- sqrt(5) * d / rho
    Phi <- (1 + scaled_d + (5 / 3) * (d / rho)^2) * exp(-scaled_d)
    
  } else {
    # General nu via modified Bessel function of the second kind (besselK)
    Phi <- matrix(1, nrow = nrow(d), ncol = ncol(d))
    nonzero <- d > 0
    scaled_d <- sqrt(2 * nu) * d[nonzero] / rho
    
    # Matérn correlation formula:
    # 2^(1 - nu) / Gamma(nu) * (sqrt(2*nu) * d / rho)^nu * K_nu(sqrt(2*nu) * d / rho)
    const <- (2^(1 - nu)) / gamma(nu)
    Phi[nonzero] <- const * (scaled_d^nu) * besselK(scaled_d, nu = nu)
  }
  
  return(Phi)
}