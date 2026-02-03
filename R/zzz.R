
# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

  requireNamespace("nimble")

  packageStartupMessage("Loading nimbleJSextras...\nRegistering multiple Jolly-Seber related distributions.\n")

  # Register the distributions explicitly for two reasons:
  # 1. Avoid message to user about automatic registrations upon first use in a nimbleModel
  # 2. Establish default len = 0 via reparameterization mechanism.

  suppressMessages({

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
        pqAvail = FALSE)), verbose = F
    )

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
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dJS_binom = list(
        BUGSdist = " dJS_binom(init, prob, size, probTrans, pstar, weight, len)",
        Rdist = " dJS_binom(init, prob, size, probTrans, pstar, weight, len)",
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
        pqAvail = FALSE)), verbose = F
    )

  })

}
