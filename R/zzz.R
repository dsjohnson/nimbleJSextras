
# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

  packageStartupMessage("Loading nimbleJSextras. Registering multiple Jolly-Seber related distributions.\n")

  # Register the distributions explicitly for two reasons:
  # 1. Avoid message to user about automatic registrations upon first use in a nimbleModel
  # 2. Establish default len = 0 via reparameterization mechanism.

  suppressMessages({

    registerDistributions(list(
      dJS_cat = list(
        BUGSdist = "dJS_cat(init, probObs, probTrans, len, pstar, checkRowSums)",
        Rdist = "dJS_cat(init, probObs, probTrans, len, pstar, checkRowSums)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'probObs = double(3)',
                  'probTrans = double(3)',
                  'len = integer(0)',
                  'pstar = double(0)',
                  'checkRowSums = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dDHMMo_bern = list(
        BUGSdist = "dDHMMo_bern(init, prob, probTrans, len)",
        Rdist = "dDHMMo_bern(init, prob, probTrans, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'prob = double(2)',
                  'probTrans = double(3)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dDHMMo_binom = list(
        BUGSdist = "dDHMMo_binom(init, prob, size, probTrans, len)",
        Rdist = "dDHMMo_binom(init, prob, size, probTrans, len)",
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
      dJS_binom = list(
        BUGSdist = " dJS_binom(init, prob, size, probTrans, len, pstar)",
        Rdist = " dJS_binom(init, prob, size, probTrans, len, pstar)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'prob = double(2)',
                  'size = double(1)',
                  'probTrans = double(3)',
                  'len = integer(0)',
                  'pstar = double(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dJS = list(
        BUGSdist = " dJS(init, prob, probTrans, len, pstar)",
        Rdist = " dJS(init, prob, probTrans, len, pstar)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'prob = double(2)',
                  'probTrans = double(3)',
                  'len = integer(0)',
                  'pstar = double(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )


  })

}
