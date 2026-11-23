## `InitMulti()` finds each replicate's files in `outputsDF`. These tests build `outputsDF` for
## two replicates that live outside `outputPath(sim)` (as they do when the run is restored from
## archives and summarised under `.../_all`) and check what InitMulti() derives from it.

analysisTimes <- c(2720, 2820, 2920, 3020)

makeOutputsDF <- function(root, reps = 1:2, repeatLastYear = 3L) {
  perRep <- function(r) {
    d <- file.path(root, paste0("rep", r))
    yrs <- c(2020, analysisTimes)
    f <- c(
      file.path(d, sprintf("vegTypeMap_year%04d.tif", yrs)),
      file.path(d, sprintf("standAgeMap_year%04d.tif", yrs)),
      file.path(d, sprintf("cohortData_year%04d.qs2", analysisTimes)),
      file.path(d, sprintf("pixelGroupMap_year%04d.tif", analysisTimes)),
      file.path(d, "rstTimeSinceFire_year2020.tif"),
      file.path(d, "flammableMap_year3020.tif")
    )
    ## the last year is registered more than once: it is saved by the sim's own save event and again
    ## by the summary modules' end-of-run save.
    c(f, rep(grep("year3020", f, value = TRUE), repeatLastYear - 1L))
  }
  ## in the order the files were saved, which is not year order for every kind of file
  set.seed(42)
  data.table::data.table(file = sample(unlist(lapply(reps, perRep))), saveTime = 3020)
}

makeMultiSim <- function(outputsDF, reps = 1:2) {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  SpaDES.core::simInit(
    times = list(start = 2020, end = 3020),
    params = list(NRV_summary = list(mode = "multi", reps = reps, summaryPeriod = c(2720L, 3020L),
                                     summaryInterval = 100L, timeSeriesTimes = 2721:2722)),
    modules = moduleName,
    objects = list(outputsDF = outputsDF, reportingPolygons = list(ELF = "placeholder")),
    paths = list(modulePath = modulePath, outputPath = file.path(root, "_all"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  )
}

runInit <- function(sim) {
  suppressWarnings( ## `LandTypeCC_reporting` is absent: asserted separately below
    SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)
  )
}

modOf <- function(sim) sim@.xData$.modObjs[[moduleName]]

test_that("each replicate contributes each analysis year once, however often a file was registered", {
  root <- withr::local_tempdir()
  m <- modOf(runInit(makeMultiSim(makeOutputsDF(root), reps = 1:2)))

  for (x in c("vtm", "sam", "cd", "pgm")) {
    expect_false(anyDuplicated(m[[x]]) > 0, label = paste("duplicated files in mod$", x))
    expect_length(m[[x]], 2L * length(analysisTimes))
  }
  ## n_reps in the envelopes is the row count per time: one file per rep and year
  years <- sub(".*_year(\\d+)\\..*", "\\1", m$vtm)
  expect_true(all(table(years) == 2L))
  expect_identical(unname(lengths(.filesByRep(m$vtm))), c(4L, 4L))
})

test_that("an analysis year is matched in the file name, not elsewhere in the path", {
  ## a temporary folder named e.g. .../file12a2820b... held "2820", so every file of every year matched
  root <- withr::local_tempdir(pattern = "run2820_")
  m <- modOf(runInit(makeMultiSim(makeOutputsDF(root), reps = 1:2)))
  for (x in c("vtm", "sam", "cd", "pgm")) {
    expect_length(m[[x]], 2L * length(analysisTimes))
  }
  expect_false(any(grepl("year2020", c(m$vtm, m$sam))))
})

test_that("the veg-type and stand-age maps of a rep are listed in the same year order", {
  ## the patch metrics pair them by position
  root <- withr::local_tempdir()
  m <- modOf(runInit(makeMultiSim(makeOutputsDF(root), reps = 1:2)))
  yearsByRep <- function(files) lapply(.filesByRep(files), function(f) sub(".*_year(\\d+)\\..*", "\\1", f))
  expect_identical(yearsByRep(m$vtm), yearsByRep(m$sam))
})

test_that("year-0 and flammable rasters come from where the replicate's outputs are", {
  root <- withr::local_tempdir()
  m <- modOf(runInit(makeMultiSim(makeOutputsDF(root), reps = 1:2)))
  rep1 <- file.path(root, "rep1")

  expect_identical(m$fvtm0, file.path(rep1, "vegTypeMap_year2020.tif"))
  expect_identical(m$fsam0, file.path(rep1, "standAgeMap_year2020.tif"))
  expect_identical(m$ftsf0, file.path(rep1, "rstTimeSinceFire_year2020.tif"))
  expect_identical(m$flm, file.path(rep1, "flammableMap_year3020.tif"))
  ## and not from outputPath(), which ends in "_all" in multi mode
  expect_false(any(grepl("_all", c(m$fvtm0, m$fsam0, m$ftsf0, m$flm))))
})

test_that("a missing LandTypeCC_reporting is warned about, not silent", {
  root <- withr::local_tempdir()
  expect_warning(
    SpaDES.core::spades(makeMultiSim(makeOutputsDF(root)), events = list(NRV_summary = "init"),
                        debug = FALSE),
    "LandTypeCC_reporting"
  )
})

test_that("national layers join the reporting polygons; standOnlyLayers get only the stand-level measures", {
  root <- withr::local_tempdir()
  sim <- makeMultiSim(makeOutputsDF(root), reps = 1:2)
  SpaDES.core::params(sim)$NRV_summary$nationalReportingLayers <- c("Province", "Protected")
  sim$studyAreaReporting <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 10, ymax = 10),
                                                                           crs = 3978)))
  asked <- NULL
  testthat::local_mocked_bindings(
    buildReportingPolygons = function(studyArea, destinationPath, targetCRS, layers, ...) {
      asked <<- layers$NAME_SHORT
      sapply(layers$NAME_SHORT, function(n) terra::vect(sf::st_sf(Name = n, geometry = sf::st_geometry(sim$studyAreaReporting))),
             simplify = FALSE)
    },
    .package = "LandWebUtils"
  )
  sim <- runInit(sim)
  m <- modOf(sim)
  expect_identical(asked, c("Province", "Protected"))
  expect_identical(names(m$reportingPolygons), c("ELF", "Province", "Protected"))
  expect_s3_class(m$reportingPolygons$Protected, "sf")
  expect_identical(names(sim$reportingPolygons), "ELF") ## the input is not changed
  expect_identical(m$standPolyNames, c("ELF", "Province", "Protected"))
  expect_identical(m$rptPolyNames, c("ELF", "Province")) ## no landscape or patch metrics of the protected areas
})

test_that("an unknown national layer stops the run before any work", {
  root <- withr::local_tempdir()
  sim <- makeMultiSim(makeOutputsDF(root), reps = 1:2)
  SpaDES.core::params(sim)$NRV_summary$nationalReportingLayers <- "Parks"
  sim$studyAreaReporting <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 10, ymax = 10),
                                                                           crs = 3978)))
  expect_error(runInit(sim), "unknown nationalReportingLayers.*Parks")
})
