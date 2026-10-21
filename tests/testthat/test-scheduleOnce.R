## `init` in mode "single" must queue map_generators and save_single once per distinct time: they
## are wanted at the start, at the start of the summary period, at the end and at each
## timeSeriesTimes year, and those times can coincide.

moduleName <- "NRV_summary"
modulePath <- dirname(normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE))
root <- withr::local_tempdir()
paths <- list(cachePath = file.path(root, "cache"), inputPath = file.path(root, "inputs"),
              modulePath = modulePath, outputPath = file.path(root, "outputs"))

singleSim <- function(summaryPeriod, timeSeriesTimes, end = 10, interval = 5L) {
  r <- terra::rast(nrows = 4, ncols = 4, vals = 1L)
  suppressMessages(SpaDES.core::simInit(
    times = list(start = 0, end = end),
    modules = moduleName,
    paths = paths,
    params = stats::setNames(
      list(list(mode = "single", summaryPeriod = summaryPeriod, summaryInterval = interval,
                timeSeriesTimes = timeSeriesTimes, .plots = NA, .useCache = FALSE)),
      moduleName),
    objects = list(cohortData = data.table::data.table(pixelGroup = 1L), pixelGroupMap = r,
                   speciesLayers = r, sppColorVect = c(Abie_bal = "red"),
                   sppEquiv = data.table::data.table(LandR = "Abie_bal"), studyAreaReporting = r)
  ))
}

queuedAfterInit <- function(...) {
  sim <- suppressMessages(SpaDES.core::spades(
    singleSim(...), events = stats::setNames(list("init"), moduleName)))
  ev <- SpaDES.core::events(sim)
  ev[ev$moduleName == moduleName, c("eventTime", "eventType")]
}

timesOf <- function(ev, type) sort(ev$eventTime[ev$eventType == type])

test_that("events are queued once per distinct time when the summary period starts at the start", {
  ev <- queuedAfterInit(c(0L, 10L), timeSeriesTimes = 0)
  expect_identical(timesOf(ev, "map_generators"), c(0, 10))
  expect_identical(timesOf(ev, "save_single"), c(0, 10))
})

test_that("events are queued once per distinct time when the summary period starts at the end", {
  ev <- queuedAfterInit(c(10L, 10L), timeSeriesTimes = 10)
  expect_identical(timesOf(ev, "map_generators"), c(0, 10))
  expect_identical(timesOf(ev, "save_single"), c(0, 10))
})

test_that("distinct times are all kept", {
  ev <- queuedAfterInit(c(5L, 10L), timeSeriesTimes = c(2, 5))
  expect_identical(timesOf(ev, "map_generators"), c(0, 2, 5, 10))
  expect_identical(timesOf(ev, "save_single"), c(0, 2, 5, 10))
})
