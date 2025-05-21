
# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

  packageStartupMessage("Loading nimbleJSextras. Registering multiple variants of the following distributions:\n ",
                        "dJS and dJS_D.\n")

  # Register the distributions explicitly for two reasons:
  # 1. Avoid message to user about automatic registrations upon first use in a nimbleModel
  # 2. Establish default len = 0 via reparameterization mechanism.

  suppressMessages({
    registerDistributions(list(
      dJS = list(
        BUGSdist = "dJS(init, probObs, probTrans, len, pstar, checkRowSums)",
        Rdist = "dJS(init, probObs, probTrans, len, pstar, checkRowSums)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'probObs = double(2)',
                  'probTrans = double(2)',
                  'len = integer(0)',
                  'pstar = double(0)',
                  'checkRowSums = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dJS_o = list(
        BUGSdist = "dJS_o(init, probObs, probTrans, len, pstar, checkRowSums)",
        Rdist = "dJS_o(init, probObs, probTrans, len, pstar, checkRowSums)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'probObs = double(3)',
                  'probTrans = double(2)',
                  'len = integer(0)',
                  'pstar = double(0)',
                  'checkRowSums = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )
    registerDistributions(list(
      dJS_D = list(
        BUGSdist = "dJS_D(init, probObs, probTrans, len, pstar, checkRowSums)",
        Rdist = "dJS_D(init, probObs, probTrans, len, pstar, checkRowSums)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'probObs = double(2)',
                  'probTrans = double(3)',
                  'len = integer(0)',
                  'pstar = double(0)',
                  'checkRowSums = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dJS_Do = list(
        BUGSdist = "dJS_Do(init, probObs, probTrans, len, pstar, checkRowSums)",
        Rdist = "dJS_Do(init, probObs, probTrans, len, pstar, checkRowSums)",
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
  })

}
