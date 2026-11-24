## Single mode checks the vegetation inputs at its first summary event, not in init: a module that
## loads after NRV_summary may provide them in its own init.

stubObjects <- function() {
  pgm <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000,
                     vals = rep(c(1, 2, 3, 1, 2), 20))
  names(pgm) <- "pixelGroup"
  cd <- data.table::data.table(pixelGroup = c(1L, 1L, 2L, 2L, 3L),
                               speciesCode = c("Pice_mar", "Betu_pap", "Betu_pap", "Pice_mar", "Pice_mar"),
                               age = c(30L, 60L, 150L, 160L, 90L), B = c(1000L, 500L, 2000L, 100L, 800L))
  sa <- terra::as.polygons(terra::ext(pgm))
  terra::crs(sa) <- terra::crs(pgm)
  list(
    cohortData = cd, pixelGroupMap = pgm, studyAreaReporting = sa,
    speciesLayers = terra::rast(pgm, vals = 1),
    sppColorVect = c(Pice_mar = "black", Betu_pap = "white", Mixed = "grey"),
    sppEquiv = data.table::data.table(LandR = c("Pice_mar", "Betu_pap"), Type = c("Conifer", "Deciduous"),
                                      EN_generic_full = c("Black spruce", "Paper birch"))
  )
}

stubSim <- function(objs, root) {
  SpaDES.core::simInit(
    times = list(start = 0, end = 2),
    params = list(NRV_summary = list(mode = "single", summaryPeriod = c(0, 2), summaryInterval = 2L,
                                timeSeriesTimes = c(0, 2), .plots = NA)),
    modules = moduleName, objects = objs,
    paths = list(modulePath = modulePath, outputPath = file.path(root, "out"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  )
}

test_that("init does not stop on inputs a later module provides", {
  objs <- stubObjects()
  objs$cohortData <- NULL
  sim <- stubSim(objs, withr::local_tempdir())
  expect_no_error(
    suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE))
  )
})

test_that("the first map_generators stops on a NULL speciesLayers", {
  objs <- stubObjects()
  objs$speciesLayers <- NULL
  sim <- stubSim(objs, withr::local_tempdir())
  expect_error(
    suppressWarnings(
      SpaDES.core::spades(sim, events = list(NRV_summary = c("init", "map_generators")), debug = FALSE)
    ),
    "needs speciesLayers in the simList"
  )
})

test_that("zero-row sppEquiv at init schedules no NRV_summary events", {
  objs <- stubObjects()
  objs$sppEquiv <- objs$sppEquiv[0L]
  sim <- stubSim(objs, withr::local_tempdir())
  expect_message(
    sim <- suppressWarnings(
      SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)
    ),
    "no tree species"
  )
  ev <- SpaDES.core::events(sim)
  expect_false(any(ev$moduleName == "NRV_summary"))
})

test_that("zero-row sppEquiv by the first event: a message, no error, nothing scheduled", {
  objs <- stubObjects()
  objs$sppEquiv <- NULL
  root <- withr::local_tempdir()
  sim <- stubSim(objs, root)
  sim <- suppressWarnings(
    SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)
  )
  ## another module's init supplies the empty table after NRV_summary's init
  sim$sppEquiv <- stubObjects()$sppEquiv[0L]
  msgs <- character()
  sim <- withCallingHandlers(
    suppressWarnings(SpaDES.core::spades(sim, debug = FALSE)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  ## two map_generators events were queued at init (times 0 and 2): the message is printed once
  expect_equal(sum(grepl("no tree species", msgs)), 1L)
  expect_length(list.files(file.path(root, "out"), pattern = "_year"), 0L)
  ## both queued times ran, and no event rescheduled itself or queued others
  done <- SpaDES.core::completed(sim)
  done <- done[done$moduleName == "NRV_summary" & done$eventType %in% c("map_generators", "save_single"), ]
  expect_equal(done$eventType, c("map_generators", "save_single", "map_generators", "save_single"))
  expect_equal(as.numeric(done$eventTime), c(0, 0, 2, 2))
})

test_that("zero-row sppEquiv with a NULL speciesLayers at the first event: a message, no error", {
  objs <- stubObjects()
  objs$sppEquiv <- NULL
  objs$speciesLayers <- NULL
  sim <- stubSim(objs, withr::local_tempdir())
  sim <- suppressWarnings(
    SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)
  )
  sim$sppEquiv <- stubObjects()$sppEquiv[0L]
  ## the no-tree check runs before the missing-input stop
  expect_message(
    expect_no_error(
      suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "map_generators"), debug = FALSE))
    ),
    "no tree species"
  )
})

test_that("an all-NA sppEquivCol is not read as 'no tree species'", {
  objs <- stubObjects()
  objs$sppEquiv$LandR <- NA_character_
  sim <- stubSim(objs, withr::local_tempdir())
  msgs <- character()
  sim <- withCallingHandlers(
    suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("no tree species", msgs)))
  expect_true(any(SpaDES.core::events(sim)$moduleName == "NRV_summary"))
})
