## vegLeadingProportion takes its default from LandR, as Biomass_core does: paramCheckOtherMods()
## stops a run whose modules disagree on it, and until 2.0.0.9028 this module's hard-coded 0.8
## disagreed with their 0.75.

paramDefault <- function(name) {
  p <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)$parameters
  p$default[[which(p$paramName == name)]]
}

test_that("vegLeadingProportion defaults to LandR::leadingSpeciesProp()", {
  expect_equal(paramDefault("vegLeadingProportion"), LandR::leadingSpeciesProp())
})

test_that("vegLeadingProportion follows the LandR.leadingSpeciesProp option", {
  withr::local_options(LandR.leadingSpeciesProp = 0.6)
  expect_equal(paramDefault("vegLeadingProportion"), 0.6)
})

## timeSeriesTimes defaults to start + 601:650, which lies past the end of a run shorter than
## 650 years. A short run gets the last 50 years of the run instead.
timeSeriesDefault <- function(end, start = 0) {
  sim <- SpaDES.core::simInit(
    times = list(start = start, end = end),
    modules = moduleName,
    paths = testPaths
  )
  SpaDES.core::params(sim)[[moduleName]]$timeSeriesTimes
}

test_that("timeSeriesTimes defaults to years inside a short run", {
  tst <- timeSeriesDefault(300)
  expect_true(all(tst >= 0 & tst <= 300))
  expect_equal(tst, 251:300)
})

test_that("timeSeriesTimes defaults to 601:650 for a 1000-year run", {
  expect_equal(timeSeriesDefault(1000), 601:650)
})

test_that("the timeSeriesTimes default clamps to a start other than 0", {
  expect_equal(timeSeriesDefault(2300, start = 2000), 2251:2300)
})

test_that("a run shorter than 50 years gets every year from start to end", {
  expect_equal(timeSeriesDefault(20), 0:20)
})

test_that("timeSeriesTimes = NA resolves to the default", {
  resolved <- function(end) {
    sim <- SpaDES.core::simInit(
      times = list(start = 0, end = end),
      params = list(NRV_summary = list(mode = "single", summaryPeriod = c(0, end),
                                       timeSeriesTimes = NA, .plots = NA)),
      modules = moduleName,
      objects = list(sppEquiv = data.table::data.table(LandR = "Pice_mar")),
      paths = testPaths
    )
    sim <- suppressWarnings(
      SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)
    )
    SpaDES.core::params(sim)[[moduleName]]$timeSeriesTimes
  }
  expect_equal(resolved(1000), 601:650)
  expect_equal(resolved(300), 251:300)
})
