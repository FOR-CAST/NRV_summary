## The mode = "multi" outputs (parquet aggregates, envelope CSVs, animation, figures) are written
## outside the simList, so they reach a targets `_files` manifest only if registered with
## registerOutputs(). Until 2.0.0.9025 none were: the stage's `_files` target tracked 0 files.

test_that("postprocess registers the files the postprocess events wrote", {
  withr::local_package("SpaDES.core") ## as at run time: module code calls SpaDES.core unqualified
  mod <- new.env()
  f <- file.path(testPaths$outputPath, c("a.parquet", "b.csv"))
  file.create(f)
  mod$ppFiles <- c(f, f[1], file.path(testPaths$outputPath, "never-written.csv"))

  register <- moduleFn(".registerPostprocessOutputs", mod)
  sim <- register(SpaDES.core::simInit(paths = testPaths))

  registered <- normalizePath(SpaDES.core::outputs(sim)$file, winslash = "/", mustWork = FALSE)
  expect_setequal(registered, normalizePath(f, winslash = "/"))
  expect_length(mod$ppFiles, 0L) ## reset, so a second call does not register twice
})

test_that("the envelope CSVs record the files they write", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  mod$ppFiles <- character(0)
  writeCSVs <- moduleFn(".writeNrvSummaryCSVs", mod)
  sim <- SpaDES.core::simInit(paths = testPaths)
  env <- data.frame(metric = c("m1", "m2"), value = c(1, 2))

  writeCSVs(sim, env, "lm_ANSR", "ANSR")

  expect_length(mod$ppFiles, 3L) ## the combined CSV + one per metric
  expect_true(all(file.exists(mod$ppFiles)))
})

test_that("the plot event keeps the sim plotFun() returns", {
  ## plotFun() registers the PNGs on the sim it returns; calling it without assigning the result
  ## discarded every figure registration.
  doEvent <- get0(paste0("doEvent.", moduleName), mode = "function")
  body <- if (is.null(doEvent)) {
    readLines(file.path(moduleRoot, paste0(moduleName, ".R")))
  } else {
    deparse(doEvent)
  }
  expect_true(any(grepl("sim <- plotFun(sim)", body, fixed = TRUE)))
})

test_that("the masked current-condition maps keep the <rep>/<prefix>_year<YYYY>.tif path nrvtools parses", {
  ## nrvtools labels each per-map result `<rep>.<prefix>_year<YYYY>_<polyName>` from the file path
  ## and stops on anything else; `_cc/cc_vegTypeMap.tif` stopped the summaries in 2.0.0.9025.
  src <- file.path(testPaths$outputPath, "maskCC", "rep01", "vegTypeMap_year0000.tif")
  dir.create(dirname(src), recursive = TRUE, showWarnings = FALSE)
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857",
                   vals = 1:4)
  terra::writeRaster(r, src, overwrite = TRUE, datatype = "INT2U")
  urbanMask <- terra::rast(r, vals = c(1L, NA, 1L, 1L)) ## cell 2 is urban

  maskCC <- moduleFn(".maskCC", new.env())
  ccDir <- file.path(testPaths$outputPath, "maskCC", "_cc")
  out <- maskCC(src, urbanMask, ccDir)

  expect_identical(normalizePath(out), normalizePath(file.path(ccDir, "rep01", "vegTypeMap_year0000.tif")))
  label <- paste0(basename(dirname(out)), ".", tools::file_path_sans_ext(basename(out)), "_ANC")
  expect_match(label, "^(.*?)\\.(.*)_year([0-9]+)_(.*)$") ## nrvtools' .parse_metric_labels() pattern
  expect_identical(terra::values(terra::rast(out), mat = FALSE), c(1L, NA, 3L, 4L))

  ## a grid mismatch leaves the map uncorrected, loudly, and returns the original path
  other <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3, crs = "EPSG:3857",
                       vals = 1L)
  expect_warning(same <- maskCC(src, other, ccDir), "does not match")
  expect_identical(same, src)
})

## makeSeralStageMapsBC() on a first replicate with the year-0 and year-1000 files, `LandTypeCC_reporting` `ccrep`
## and the parameter `sieveThresh`. The
## event needs BC's inputs (ecoregionLayer polygons with Natural Disturbance Types, BC species) and a callr plan, so
## nrvtools::writeSeralStageMapBC() is replaced by a stub that writes a seral map of four classes on the grid of each
## pixel-group map, and P() and plan() by stubs too. Returns the module state and the root folder.
bcSeralRun <- function(ccrep, sieveThresh = 2L, .env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .env)
  out <- file.path(root, "mainSim")
  rep1 <- file.path(out, "rep01")
  dir.create(rep1, recursive = TRUE)
  crs <- "EPSG:3400"
  rtm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = crs)
  for (y in c("0000", "1000")) {
    terra::writeRaster(terra::rast(rtm, vals = 1:4), file.path(rep1, paste0("pixelGroupMap_year", y, ".tif")))
    file.create(file.path(rep1, paste0("cohortData_year", y, ".qs2")))
  }
  poly <- terra::as.polygons(terra::ext(rtm), crs = crs)
  ndt <- sf::st_as_sf(poly)
  ndt$NDTBEC <- "NDT3_SBS"
  sim <- SpaDES.core::simInit(
    paths = list(outputPath = out, cachePath = file.path(root, "cache"), inputPath = file.path(root, "in")),
    objects = c(list(studyAreaReporting = poly, reportingPolygons = list(ecoregionLayer = ndt)),
                if (!is.null(ccrep)) list(LandTypeCC_reporting = ccrep))
  )
  mod <- new.env()
  mod$rep1Dir <- c("1" = rep1)
  mod$cohortFiles <- file.path(rep1, c("cohortData_year0000.qs2", "cohortData_year1000.qs2"))
  mod$cd <- mod$cohortFiles[2]
  mod$pgm <- file.path(rep1, "pixelGroupMap_year1000.tif")
  bc <- moduleFn("makeSeralStageMapsBC", mod)
  e <- environment(bc)
  e$P <- function(sim, ...) list(simTimes = c(0, 1000), sieveThresh = sieveThresh)
  e$plan <- function(...) NULL
  e$rast <- terra::rast ## attached in the module's environment
  e$writeSeralStageMapBC <- function(cd, pgm, ndtbec, ssm = sub("pixelGroupMap", "seralStageMap", pgm)) {
    for (i in seq_along(pgm)) {
      r <- terra::rast(terra::rast(pgm[i]), vals = c(1L, 1L, 2L, 2L)) ## two clumps of 2, which the sieve keeps
      levels(r) <- data.frame(id = 1:4, values = c("early", "mid", "mature", "old"))
      terra::writeRaster(r, ssm[i], datatype = "INT1U", overwrite = TRUE)
    }
    ssm
  }
  bc(sim)
  list(mod = mod, rep1 = rep1, out = out, rtm = rtm)
}

test_that("bc's current-condition seral map leaves out the urban pixels of LandTypeCC_reporting", {
  withr::local_package("SpaDES.core") ## as at run time: module code calls SpaDES.core unqualified
  rtm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = "EPSG:3400")
  ccrep <- terra::rast(rtm, vals = c(1L, 17L, 17L, 5L)) ## cells 2 and 3 urban (17, LCC 2020)
  r <- bcSeralRun(ccrep)
  ssm0 <- unname(r$mod$ssm0)
  expect_identical(ssm0, file.path(r$out, "_cc", "rep01", "seralStageMap-sieve2_year0000.tif"))
  expect_equal(terra::values(terra::rast(ssm0), mat = FALSE), c(1, NA, NA, 2))
  ## the replicates' maps are not masked: the envelope has no urban land
  expect_identical(unname(r$mod$ssm), file.path(r$rep1, "seralStageMap-sieve2_year1000.tif"))
  expect_false(anyNA(terra::values(terra::rast(r$mod$ssm), mat = FALSE)))

  ## without LandTypeCC_reporting (InitMulti() has warned): the map is not masked
  r0 <- bcSeralRun(NULL)
  expect_identical(unname(r0$mod$ssm0), file.path(r0$rep1, "seralStageMap-sieve2_year0000.tif"))
  expect_false(anyNA(terra::values(terra::rast(r0$mod$ssm0), mat = FALSE)))
})

test_that("bc sieves only with a sieveThresh of 2 or more, and its default does not", {
  ## terra::sieve() stops on a threshold below 2; a threshold of 1 merges nothing
  withr::local_package("SpaDES.core")
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  default <- md$parameters$default[[which(md$parameters$paramName == "sieveThresh")]]
  unsieved <- function(r) file.path(r$rep1, c("seralStageMap_year0000.tif", "seralStageMap_year1000.tif"))
  for (thr in list(default, 1L, NA_integer_)) {
    r <- bcSeralRun(NULL, sieveThresh = thr)
    expect_identical(unname(c(r$mod$ssm0, r$mod$ssm)), unsieved(r), info = paste("sieveThresh", thr))
    expect_identical(list.files(r$rep1, "sieve"), character(0), info = paste("sieveThresh", thr))
  }
  r2 <- bcSeralRun(NULL, sieveThresh = 2L)
  sieved <- file.path(r2$rep1, c("seralStageMap-sieve2_year0000.tif", "seralStageMap-sieve2_year1000.tif"))
  expect_identical(unname(c(r2$mod$ssm0, r2$mod$ssm)), sieved)
  expect_true(all(file.exists(sieved)))
})

test_that("bc's sieved seral maps get through the patch metrics, which take the year from their names", {
  ## nrvtools labels each map's result `<rep>.<file name without .tif>_<poly>` and stops unless a `_` follows the
  ## year's digits; `seralStageMap_year0000-sieve2` stopped patchMetricsSeralBC() ("cannot parse metric label(s)")
  withr::local_package("SpaDES.core")
  withr::local_package("nrvtools") ## as at run time: the event calls nrvtools unqualified
  rtm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = "EPSG:3400")
  ccrep <- terra::rast(rtm, vals = c(1L, 17L, 17L, 5L)) ## cells 2 and 3 urban: the masked map is in `_cc/rep01/`
  r <- bcSeralRun(ccrep, sieveThresh = 2L)

  ## the event's patch metrics on those maps, for a reporting layer named NDTBEC (which also writes SeralTable.csv)
  ndt <- sf::st_as_sf(terra::as.polygons(terra::ext(rtm), crs = "EPSG:3400"))
  ndt$NDTBEC <- "NDT3_SBS"
  attr(ndt, "field") <- "NDTBEC"
  root <- dirname(r$out)
  sim <- SpaDES.core::simInit(
    paths = list(outputPath = r$out, cachePath = file.path(root, "cache"), inputPath = file.path(root, "in")),
    objects = list(studyAreaReporting = terra::as.polygons(terra::ext(rtm), crs = "EPSG:3400"),
                   reportingPolygons = list(NDTBEC = ndt))
  )
  r$mod$flm <- file.path(root, "flammableMap.tif")
  terra::writeRaster(terra::rast(rtm, vals = 1L), r$mod$flm)
  pm <- moduleFn("patchMetricsSeralBC", r$mod)
  e <- environment(pm)
  e$Cache <- function(FUN, ..., .cacheExtra = NULL) FUN(...)
  build <- moduleFn(".buildRepDataset", r$mod)
  environment(build)$.summarizeRepDataset <- moduleFn(".summarizeRepDataset", r$mod)
  e$.buildRepDataset <- build
  e$.writeNrvSummaryCSVs <- moduleFn(".writeNrvSummaryCSVs", r$mod)

  expect_no_error(suppressMessages(pm(sim)))
  expect_equal(unique(r$mod$sspm_NDTBEC_CC$time), 0)
  expect_equal(unique(r$mod$sspm_NDTBEC$time), 1000)
  seral <- utils::read.csv(file.path(r$out, "SeralTable.csv"))
  expect_setequal(seral$class, c("early", "mid"))

  ## nrvtools' own parser on the labels it builds from the maps' paths
  ssm <- unname(c(r$mod$ssm0, r$mod$ssm))
  lbl <- paste0(basename(dirname(ssm)), ".", tools::file_path_sans_ext(basename(ssm)), "_NDT3_SBS")
  expect_identical(nrvtools:::.parse_metric_labels(lbl),
                   list(rep = c(1L, 1L), time = c(0L, 1000L), poly = c("NDT3_SBS", "NDT3_SBS")))
})
