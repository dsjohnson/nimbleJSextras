
# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

  packageStartupMessage("Loading nimbleJSextras...\nRegistering multiple Jolly-Seber related distributions.\n")

  # Register the distributions explicitly for two reasons:
  # 1. Avoid message to user about automatic registrations upon first use in a nimbleModel
  # 2. Establish default len = 0 via reparameterization mechanism.

  suppressMessages({

    registerDistributions(list(
      dJS = list(
        BUGSdist = " dJS(pstar, piVector, PArray, GammaArray, len)",
        Rdist = " dJS(pstar, piVector, PArray, GammaArray, len=0)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'pstar = double(0)',
                  'piVector = double(1)',
                  'PArray = double(2)',
                  'GammaArray = double(3)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dJS_ms = list(
        BUGSdist = "dJS_ms(pstar, piVector, PArray, GammaArray, len)",
        Rdist = "dJS_ms(pstar, piVector, PArray, GammaArray, len=0)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'pstar = double(0)',
                  'piVector= double(1)',
                  'PArray = double(3)',
                  'GammaArray = double(3)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dJS_rd = list(
        BUGSdist = " dJS_rd(pstar, piVector, PArray, nSubOcc, GammaArray, len)",
        Rdist = " dJS_rd(pstar, piVector, PArray, nSubOcc, GammaArray, len=0)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'pstar = double(0)',
                  'piVector = double(1)',
                  'PArray = double(2)',
                  'nSubOcc = double(1)',
                  'GammaArray = double(3)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dhmm_bern = list(
        BUGSdist = "dhmm_bern(piVector, PArray, GammaArray, len)",
        Rdist = "dhmm_bern(piVector, PArray, GammaArray, len=0)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'piVector= double(1)',
                  'PArray = double(3)',
                  'GammaArray = double(3)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dhmm_binom = list(
        BUGSdist = "dhmm_binom(piVector, PArray, size, GammaArray, len)",
        Rdist = "dhmm_binom(piVector, PArray, size, GammaArray, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'piVector= double(1)',
                  'PArray = double(3)',
                  'size = double(1)',
                  'GammaArray = double(3)',
                  'len = integer(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dState = list(
        BUGSdist = " dState(piVector, PArray, GammaArray, x)",
        Rdist = " dState(piVector, PArray, GammaArray, x)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'piVector = double(1)',
                  'PArray = double(2)',
                  'GammaArray = double(3)',
                  'x_cond = double(1)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dState_ms = list(
        BUGSdist = " dState_ms(piVector, PArray, GammaArray, x)",
        Rdist = " dState_ms(piVector, PArray, GammaArray, x)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'piVector = double(1)',
                  'PArray = double(2)',
                  'GammaArray = double(3)',
                  'x_cond = double(1)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dState_rd = list(
        BUGSdist = " dState_rd(piVector, PArray, nSubOcc, GammaArray, x)",
        Rdist = " dState_rd(piVector, PArray, nSubOcc, GammaArray, x)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'piVector = double(1)',
                  'PArray = double(2)',
                  'nSubOcc = double(1)',
                  'GammaArray = double(3)',
                  'x_cond = double(1)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )

    registerDistributions(list(
      dnu = list(
        BUGSdist = " dnu(nUndet, piVector, PArray, GammaArray, len)",
        Rdist = " dnu(nUndet, piVector, PArray, GammaArray, len)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'nUndet = integer(0)',
                  'piVector = double(1)',
                  'PArray = double(2)',
                  'GammaArray = double(3)',
                  'len = double(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE)), verbose = F
    )



  })

}
