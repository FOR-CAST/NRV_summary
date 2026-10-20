## Year 0 is the current condition that mode = "multi" reports today's values from. Through
## 2.0.0.9029 it was saved at .last(), after every other time-0 event: with Biomass_core's
## `growthInitialTime = start(sim)` (LandWeb's setting), it held a simulated year of growth and
## old-age mortality, so today's values moved with the species traits.

## A stand-in for Biomass_core's time-0 growth: at start(sim), at Biomass_core's priority for
## mortalityAndGrowth (6), it doubles every cohort's biomass.
writeGrowthStandIn <- function(dir) {
  d <- file.path(dir, "growthStandIn")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    'defineModule(sim, list(',
    '  name = "growthStandIn", description = "test stand-in for a time-0 growth event",',
    '  keywords = "test", authors = person("Test", "Author"), childModules = character(0),',
    '  version = list(growthStandIn = "0.0.1"), timeframe = as.POSIXlt(c(NA, NA)),',
    '  timeunit = "year", citation = list(), documentation = list(), reqdPkgs = list("data.table"),',
    '  parameters = rbind(),',
    '  inputObjects = bindrows(expectsInput("cohortData", "data.table", "cohorts")),',
    '  outputObjects = bindrows(createsOutput("cohortData", "data.table", "cohorts, grown"))',
    '))',
    'doEvent.growthStandIn <- function(sim, eventTime, eventType) {',
    '  switch(eventType,',
    '    init = sim <- scheduleEvent(sim, start(sim), "growthStandIn", "grow", 6),',
    '    grow = {',
    '      cd <- data.table::copy(sim$cohortData)',
    '      data.table::set(cd, j = "B", value = cd$B * 2L)',
    '      sim$cohortData <- cd',
    '    }',
    '  )',
    '  invisible(sim)',
    '}'
  ), file.path(d, "growthStandIn.R"))
  dir
}

test_that("year 0 is saved before any time-0 dynamics, once", {
  mods <- writeGrowthStandIn(withr::local_tempdir())
  out <- file.path(testPaths$outputPath, "currentCondition")
  crs <- "EPSG:3400"
  rtm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = crs)
  spp <- c("Pice_gla", "Popu_tre")
  sppEquiv <- LandR::sppEquivalencies_CA
  sppEquiv <- sppEquiv[sppEquiv[["LandR"]] %in% spp, ]
  cohortData <- data.table::data.table(
    pixelGroup = c(1L, 1L, 2L, 3L),
    speciesCode = factor(c("Pice_gla", "Popu_tre", "Pice_gla", "Popu_tre"), levels = spp),
    age = c(80L, 40L, 150L, 20L),
    B = c(6000L, 2000L, 9000L, 1500L)
  )
  speciesLayers <- c(terra::rast(rtm, vals = c(60L, 60L, 100L, 0L)),
                     terra::rast(rtm, vals = c(40L, 40L, 0L, 100L)))
  names(speciesLayers) <- spp

  sim <- SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 10),
    params = list(NRV_summary = list(
      mode = "single", summaryPeriod = c(5L, 10L), summaryInterval = 5L,
      timeSeriesTimes = c(0, 1), ## year 0 among them must not add a second, later year-0 save
      ccStandAgeSource = "cohortData" ## no age layer here; see the stand-age tests below
    )),
    modules = c("NRV_summary", "growthStandIn"),
    objects = list(
      cohortData = data.table::copy(cohortData),
      pixelGroupMap = terra::rast(rtm, vals = c(1L, 1L, 2L, 3L)),
      speciesLayers = speciesLayers,
      sppEquiv = sppEquiv,
      sppColorVect = LandR::sppColors(sppEquiv, "LandR", newVals = "Mixed", palette = "Accent"),
      studyAreaReporting = terra::as.polygons(terra::ext(rtm), crs = crs),
      flammableMap = terra::rast(rtm, vals = 1L)
    ),
    paths = utils::modifyList(testPaths, list(modulePath = c(testPaths$modulePath, mods),
                                              outputPath = out))
  )

  ## the stand-in ran at time 0, after the year-0 save
  expect_identical(sim$cohortData$B, cohortData$B * 2L)
  saved <- qs2::qs_read(file.path(out, "cohortData_year00.qs2"))
  expect_identical(saved$B, cohortData$B)
  cmp <- as.data.frame(SpaDES.core::completed(sim))
  expect_identical(sum(cmp$eventType == "save_single" & as.numeric(cmp$eventTime) == 0), 1L)
})

## The current condition's stand ages. The cohorts' ages are binned when they are made, and
## LandR::standAgeMapGenerator() floors each pixel's mean age to the decade, so through 2.0.0.9033
## the year-0 map moved part of today's forest into older age classes than the age layer the
## cohorts were made from.

## A 2 x 2 landscape and an age layer one cell larger on each side, aligned with it. Cell 2 has
## cohorts but no value in the age layer, cell 4 an age but no cohorts. By the cohorts (weighted by
## biomass, floored to the decade), the ages are 70, 70, 150 and NA.
ccAgeFixture <- function() {
  crs <- "EPSG:3400"
  rtm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = crs)
  ages <- rep(999, 16)
  ages[c(6, 7, 10, 11)] <- c(86.5, NA, 12.5, 300) ## the cells of `rtm`; round() gives 86 and 12
  spp <- c("Pice_gla", "Popu_tre")
  sppEquiv <- LandR::sppEquivalencies_CA
  sppEquiv <- sppEquiv[sppEquiv[["LandR"]] %in% spp, ]
  speciesLayers <- c(terra::rast(rtm, vals = c(75L, 75L, 100L, 0L)),
                     terra::rast(rtm, vals = c(25L, 25L, 0L, 0L)))
  names(speciesLayers) <- spp
  list(
    cohortData = data.table::data.table(
      pixelGroup = c(1L, 1L, 2L),
      speciesCode = factor(c("Pice_gla", "Popu_tre", "Pice_gla"), levels = spp),
      age = c(80L, 40L, 150L),
      B = c(6000L, 2000L, 9000L)
    ),
    pixelGroupMap = terra::rast(rtm, vals = c(1L, 1L, 2L, NA)),
    standAgeMap = terra::rast(terra::extend(rtm, 1), vals = ages),
    speciesLayers = speciesLayers,
    sppEquiv = sppEquiv,
    sppColorVect = LandR::sppColors(sppEquiv, "LandR", newVals = "Mixed", palette = "Accent"),
    studyAreaReporting = terra::as.polygons(terra::ext(rtm), crs = crs),
    flammableMap = terra::rast(rtm, vals = 1L)
  )
}

runCurrentConditionAge <- function(objects, params, dir) {
  defaults <- list(mode = "single", summaryPeriod = c(5L, 10L), summaryInterval = 5L,
                   timeSeriesTimes = c(0, 1))
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 10),
    params = list(NRV_summary = utils::modifyList(defaults, params)),
    modules = "NRV_summary",
    objects = objects,
    paths = utils::modifyList(testPaths, list(outputPath = dir))
  )
}

savedAges <- function(dir, year) {
  f <- file.path(dir, sprintf("standAgeMap_year%02d.tif", year))
  terra::values(terra::rast(f), mat = FALSE)
}

test_that(".ccStandAgeMap() takes the age layer, rounded half up, on the cells with cohorts", {
  f <- ccAgeFixture()
  ccStandAgeMap <- moduleFn(".ccStandAgeMap", new.env())
  cohortAge <- LandR::standAgeMapGenerator(f$cohortData, f$pixelGroupMap)
  expect_equal(terra::values(cohortAge, mat = FALSE), c(70, 70, 150, NA))

  expect_message(out <- ccStandAgeMap(f$standAgeMap, cohortAge), "1 pixels with cohorts have no")
  expect_true(terra::compareGeom(out, cohortAge))
  expect_equal(terra::values(out, mat = FALSE), c(87, 70, 13, NA))
})

test_that(".ccStandAgeMap() warns and keeps the cohort ages without an aligned age layer", {
  f <- ccAgeFixture()
  ccStandAgeMap <- moduleFn(".ccStandAgeMap", new.env())
  cohortAge <- LandR::standAgeMapGenerator(f$cohortData, f$pixelGroupMap)

  expect_warning(out <- ccStandAgeMap(NULL, cohortAge), "is absent")
  expect_equal(terra::values(out, mat = FALSE), c(70, 70, 150, NA))

  shifted <- terra::shift(f$standAgeMap, dx = 100) ## cells are 250 m
  expect_warning(out <- ccStandAgeMap(shifted, cohortAge), "not on a grid aligned")
  expect_equal(terra::values(out, mat = FALSE), c(70, 70, 150, NA))
})

test_that("year 0's stand ages come from the age layer, the later years' from the cohorts", {
  f <- ccAgeFixture()
  out <- file.path(testPaths$outputPath, "currentConditionAge", "standAgeMap")
  sim <- runCurrentConditionAge(f, list(), out)

  expect_equal(savedAges(out, 0), c(87, 70, 13, NA))
  for (year in c(1, 5, 10)) {
    expect_equal(savedAges(out, year), c(70, 70, 150, NA))
  }
  ## the input is read, not changed
  expect_equal(terra::values(sim$standAgeMap, mat = FALSE),
               terra::values(f$standAgeMap, mat = FALSE))
})

test_that("ccStandAgeSource = 'cohortData' keeps the cohort ages at year 0", {
  f <- ccAgeFixture()
  out <- file.path(testPaths$outputPath, "currentConditionAge", "cohortData")
  runCurrentConditionAge(f, list(ccStandAgeSource = "cohortData"), out)

  expect_equal(savedAges(out, 0), c(70, 70, 150, NA))
})
