# #' @importFrom utils browseURL
# #' @export
# find_examples <- function(x="dipper", open=TRUE){
#   if(x=="dipper"){
#     ddd <- system.file("examples","dipper_js", package="nimbleJSextras")
#   } else if(x=="honu"){
#     ddd <- system.file("examples","honu_js", package="nimbleJSextras")
#   }
#   if(open){
#     utils::browseURL(ddd)
#   } else{
#       return(ddd)
#     }
# }

###
### R/nimble_types.R
###
#' @title Simulation output object
#' @field ch A matrix of observed states for the simulated data (included all zero capture hostories as well)
#' @field true_states A matrix of true states for the simulated data
#' @field observed An indicator vector that shows which individuals were actually captured.
#' @import nimble
#' @export
sim_js_list <- nimbleList(
  ch = double(2),
  true_states = double(2),
  observed = integer(1)
)

###
### General packages to import
###
#' @import methods
#' @importFrom stats dbinom dpois rbinom rpois

###
### This registers all distributions when the package is loaded.
###
.onAttach <- function(libname, pkgname) {

  requireNamespace("nimble")

  packageStartupMessage("Loading nimbleJSextras...\nRegistering multiple Jolly-Seber related distributions.\n")

  # Register the distributions explicitly for two reasons:
  # 1. Avoid message to user about automatic registrations upon first use in a nimbleModel
  # 2. Establish default len = 0 via reparameterization mechanism.

  suppressMessages({

    # -------------------------------------------------------------------------
    # Standard HMM Distributions
    # -------------------------------------------------------------------------
    registerDistributions(list(
      dhmm_binom = list(
        BUGSdist = "dhmm_binom(init, prob, size, probTrans, len)",
        Rdist = "dhmm_binom(init, prob, size, probTrans, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'prob = double(2)',
                  'size = double(1)',
                  'probTrans = double(3)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dhmm_pois = list(
        BUGSdist = "dhmm_pois(init, rate, probTrans, len)",
        Rdist = "dhmm_pois(init, rate, probTrans, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'rate = double(2)',
                  'probTrans = double(3)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    # -------------------------------------------------------------------------
    # Jolly-Seber Marginal Likelihood Distributions
    # -------------------------------------------------------------------------
    registerDistributions(list(
      dJS_ms = list(
        BUGSdist = "dJS_ms(init, probObs, probTrans, pstar, weight, len)",
        Rdist = "dJS_ms(init, probObs, probTrans, pstar, weight, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'probObs = double(3)',
                  'probTrans = double(3)',
                  'pstar = double(0)',
                  'weight = double(0)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJS_binom = list(
        BUGSdist = "dJS_binom(init, prob, size, probTrans, pstar, weight, len)",
        Rdist = "dJS_binom(init, prob, size, probTrans, pstar, weight, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'prob = double(2)',
                  'size = double(1)',
                  'probTrans = double(3)',
                  'pstar = double(0)',
                  'weight = double(0)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJS_pois = list(
        BUGSdist = "dJS_pois(init, rate, probTrans, pstar, weight, len)",
        Rdist = "dJS_pois(init, rate, probTrans, pstar, weight, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'rate = double(2)',
                  'probTrans = double(3)',
                  'pstar = double(0)',
                  'weight = double(0)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJS_mark = list(
        BUGSdist = "dJS_mark(mark, init, prob, probTrans, pstar, knots, sigma_knot, len)",
        Rdist = "dJS_mark(mark, init, prob, probTrans, pstar, knots, sigma_knot, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'mark = double(1)',
                  'init = double(1)',
                  'prob = double(2)',
                  'probTrans = double(3)',
                  'pstar = double(0)',
                  'knots = double(1)',
                  'sigma_knot = double()',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    # -------------------------------------------------------------------------
    # Observed Conditional State Distributions (dJSCond_*)
    # -------------------------------------------------------------------------
    registerDistributions(list(
      dJSCond_ms = list(
        BUGSdist = "dJSCond_ms(obs, init, probObs, probTrans)",
        Rdist = "dJSCond_ms(obs, init, probObs, probTrans)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'obs = double(1)',
                  'init = double(1)',
                  'probObs = double(3)',
                  'probTrans = double(3)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJSCond_binom = list(
        BUGSdist = "dJSCond_binom(obs, init, prob, size, probTrans)",
        Rdist = "dJSCond_binom(obs, init, prob, size, probTrans)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'obs = double(1)',
                  'init = double(1)',
                  'prob = double(2)',
                  'size = double(1)',
                  'probTrans = double(3)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJSCond_pois = list(
        BUGSdist = "dJSCond_pois(obs, init, rate, probTrans)",
        Rdist = "dJSCond_pois(obs, init, rate, probTrans)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'obs = double(1)',
                  'init = double(1)',
                  'rate = double(2)',
                  'probTrans = double(3)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJSCond_mark = list(
        BUGSdist = "dJSCond_mark(obs, mark, init, prob, probTrans, knots, sigma_knot)",
        Rdist = "dJSCond_mark(obs, mark, init, prob, probTrans, knots, sigma_knot)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'obs = double(1)',
                  'mark = double(1)',
                  'init = double(1)',
                  'prob = double(2)',
                  'probTrans = double(3)',
                  'knots = double(1)',
                  'sigma_knot = double()'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    # -------------------------------------------------------------------------
    # Unobserved Abundance State Distributions (dJSUndet_*)
    # -------------------------------------------------------------------------
    registerDistributions(list(
      dJSUndet_ms = list(
        BUGSdist = "dJSUndet_ms(nu, init, probObs, probTrans)",
        Rdist = "dJSUndet_ms(nu, init, probObs, probTrans)",
        discrete = TRUE,
        types = c('value = double(2)',
                  'nu = double(0)',
                  'init = double(1)',
                  'probObs = double(3)',
                  'probTrans = double(3)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJSUndet_binom = list(
        BUGSdist = "dJSUndet_binom(nu, init, prob, size, probTrans)",
        Rdist = "dJSUndet_binom(nu, init, prob, size, probTrans)",
        discrete = TRUE,
        types = c('value = double(2)',
                  'nu = double(0)',
                  'init = double(1)',
                  'prob = double(2)',
                  'size = double(1)',
                  'probTrans = double(3)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJSUndet_pois = list(
        BUGSdist = "dJSUndet_pois(nu, init, rate, probTrans)",
        Rdist = "dJSUndet_pois(nu, init, rate, probTrans)",
        discrete = TRUE,
        types = c('value = double(2)',
                  'nu = double(0)',
                  'init = double(1)',
                  'rate = double(2)',
                  'probTrans = double(3)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

    registerDistributions(list(
      dJSUndet_mark = list(
        BUGSdist = "dJSUndet_mark(nu, init, prob, probTrans)",
        Rdist = "dJSUndet_mark(nu, init, prob, probTrans)",
        discrete = TRUE,
        types = c('value = double(2)',
                  'nu = double(0)',
                  'init = double(1)',
                  'prob = double(2)',
                  'probTrans = double(3)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = FALSE
    )

  })

}
