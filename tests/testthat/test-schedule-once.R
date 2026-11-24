## In mode = "single", `init` queues map_generators and save_single at the start, at the start of
## the summary period, at the end and at each timeSeriesTimes year, and each also reschedules itself
## summaryInterval later inside the summary period. Each distinct time must run once, and year 0
## keeps its place after every module's init and before time-0 growth (.first() + 1).

test_that("each summary time has one map_generators and one save_single", {
  out <- file.path(testPaths$outputPath, "scheduleOnce")
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
    times = list(start = 0, end = 1000),
    params = list(NRV_summary = list(
      mode = "single", summaryPeriod = c(700, 1000), summaryInterval = 50,
      timeSeriesTimes = 601:650, ccStandAgeSource = "cohortData"
    )),
    modules = "NRV_summary",
    objects = list(
      cohortData = cohortData,
      pixelGroupMap = terra::rast(rtm, vals = c(1L, 1L, 2L, 3L)),
      speciesLayers = speciesLayers,
      sppEquiv = sppEquiv,
      sppColorVect = LandR::sppColors(sppEquiv, "LandR", newVals = "Mixed", palette = "Accent"),
      studyAreaReporting = terra::as.polygons(terra::ext(rtm), crs = crs),
      flammableMap = terra::rast(rtm, vals = 1L)
    ),
    paths = utils::modifyList(testPaths, list(outputPath = out))
  )

  cmp <- as.data.frame(SpaDES.core::completed(sim))
  cmp$eventTime <- as.numeric(cmp$eventTime)
  want <- c(0, 601:650, seq(700, 1000, 50))
  for (type in c("map_generators", "save_single")) {
    times <- cmp$eventTime[cmp$eventType == type]
    expect_identical(sort(unique(times)), sort(as.numeric(want)), info = type)
    expect_equal(as.vector(table(times)), rep(1L, length(want)), info = type)
  }

  y0 <- cmp[cmp$eventTime == 0 & cmp$eventType %in% c("map_generators", "save_single"), ]
  expect_equal(nrow(y0), 2L)
  expect_true(all(y0$eventPriority == SpaDES.core::.first() + 1))
})
