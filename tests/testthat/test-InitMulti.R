## InitMulti() finds each replicate's saved files in the `outputsDF` object when the run registered them
## (FireSense: folders `rep1`, `rep2`, ...) and on disk under `outputPath(sim)` when it did not
## (LandWeb: folders `rep01`, `rep02`, ..., nothing registered).

yrs <- c(700, 1000)

## the files InitMulti() expects of each replicate, as paths under `d`
repFiles <- function(d, years = yrs) {
  c(
    file.path(d, sprintf("cohortData_year%04d.qs2", years)),
    file.path(d, sprintf("pixelGroupMap_year%04d.tif", years)),
    file.path(d, sprintf("standAgeMap_year%04d.tif", years)),
    file.path(d, sprintf("vegTypeMap_year%04d.tif", years)),
    file.path(d, "flammableMap_year1000.tif")
  )
}

## write the files, empty
touchFiles <- function(files) {
  for (d in unique(dirname(files))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.create(files)
  invisible(files)
}

initMultiSim <- function(root, reps, outputsFiles = NULL) {
  sim <- SpaDES.core::simInit(
    times = list(start = 0, end = 1000),
    params = list(NRV_summary = list(mode = "multi", reps = reps, simTimes = c(0, 1000), summaryPeriod = c(700L, 1000L),
                                     summaryInterval = 300L, timeSeriesTimes = 701:702)),
    modules = moduleName,
    objects = c(list(reportingPolygons = list(ELF = "placeholder")),
                if (!is.null(outputsFiles)) list(outputsDF = data.frame(file = outputsFiles))),
    paths = list(modulePath = modulePath, outputPath = file.path(root, "mainSim"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  )
  sim
}

runInit <- function(sim) {
  suppressWarnings( ## `LandTypeCC_reporting` is absent
    SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)
  )
}

modOf <- function(sim) sim@.xData$.modObjs[[moduleName]]

test_that("rep01 folders, no outputsDF: all reps found", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, sprintf("rep%02d", 1:3)), repFiles)))

  m <- modOf(runInit(initMultiSim(root, reps = 1:3)))
  expect_identical(m$allReps, c("rep01", "rep02", "rep03"))
  expect_identical(m$rep1Dir, file.path(out, "rep01"))
  expect_identical(m$flm, file.path(out, "rep01", "flammableMap_year1000.tif"))
  expect_length(m$vtm, 3L * length(yrs))
  expect_length(m$cd, 3L * length(yrs))
})

test_that("rep1 folders listed in outputsDF", {
  root <- withr::local_tempdir()
  ## the replicates are outside outputPath(), and the last year is registered twice
  files <- unlist(lapply(file.path(root, c("rep1", "rep2")), repFiles))
  files <- c(files, grep("year1000", files, value = TRUE))
  touchFiles(files)

  m <- modOf(runInit(initMultiSim(root, reps = 1:2, outputsFiles = files)))
  expect_identical(m$allReps, c("rep1", "rep2"))
  expect_identical(m$rep1Dir, file.path(root, "rep1"))
  expect_identical(m$fvtm0, file.path(root, "rep1", "vegTypeMap_year0000.tif"))
  expect_identical(m$flm, file.path(root, "rep1", "flammableMap_year1000.tif"))
  for (x in c("vtm", "sam", "cd", "pgm")) {
    expect_false(anyDuplicated(m[[x]]) > 0, label = paste("duplicated files in mod$", x))
    expect_length(m[[x]], 2L * length(yrs))
  }
})

test_that("a missing file stops InitMulti with or without outputsDF", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  files <- unlist(lapply(file.path(out, sprintf("rep%02d", 1:3)), repFiles))
  gone <- file.path(out, "rep02", "vegTypeMap_year1000.tif")
  touchFiles(setdiff(files, gone))
  expect_error(runInit(initMultiSim(root, reps = 1:3)), "1 simulation files appear to be missing")

  root2 <- withr::local_tempdir()
  files2 <- unlist(lapply(file.path(root2, c("rep1", "rep2", "rep3")), repFiles))
  gone2 <- file.path(root2, "rep2", "vegTypeMap_year1000.tif")
  touchFiles(setdiff(files2, gone2))
  expect_error(
    runInit(initMultiSim(root2, reps = 1:3, outputsFiles = setdiff(files2, gone2))),
    "1 simulation files appear to be missing"
  )
})

test_that("reps = 1:2 ignores rep10-rep12 on disk", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, sprintf("rep%02d", c(1:2, 10:12))), repFiles)))

  m <- modOf(runInit(initMultiSim(root, reps = 1:2)))
  expect_identical(m$allReps, c("rep01", "rep02"))
  expect_false(any(grepl("rep1[0-2]", c(m$vtm, m$sam, m$cd, m$pgm))))
  expect_length(m$vtm, 2L * length(yrs))
})

test_that("fewer reps found than asked stops", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, sprintf("rep%02d", 1:2)), repFiles)))
  expect_error(runInit(initMultiSim(root, reps = 1:3)), "replicate")

  ## and with nothing at all
  root0 <- withr::local_tempdir()
  expect_error(runInit(initMultiSim(root0, reps = 1:3)), "replicate")
})

test_that("file-type filters read the file name, not the folder names", {
  ## a folder named with `png` or `txt` neither drops files, and a sidecar `.xml` is still left out
  root <- withr::local_tempdir(pattern = "png_txt_")
  out <- file.path(root, "mainSim")
  files <- unlist(lapply(file.path(out, sprintf("rep%02d", 1:2)), repFiles))
  touchFiles(c(files, file.path(out, "rep01", "vegTypeMap_year0700.tif.aux.xml")))

  m <- modOf(runInit(initMultiSim(root, reps = 1:2)))
  expect_false(any(grepl("xml$", m$vtm)))
  expect_length(m$vtm, 2L * length(yrs))
  expect_length(m$sam, 2L * length(yrs))
})

test_that(".discoverOutputFiles uses the outputsDF when it has rows, the disk otherwise", {
  discover <- moduleFn(".discoverOutputFiles", new.env())
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(file.path(out, "rep01", c("b.tif", "a.tif")))
  sim <- initMultiSim(root, reps = 1L)
  expect_setequal(basename(discover(sim)), c("a.tif", "b.tif"))

  sim$outputsDF <- data.frame(file = file.path(root, "z", c("q_year0001.tif", "p_year0001.tif", "q_year0001.tif")))
  expect_identical(discover(sim), file.path(root, "z", c("p_year0001.tif", "q_year0001.tif")))
})

test_that("BC seral and animation read the first replicate folder found", {
  rep1 <- file.path("x", "rep1")
  year0 <- moduleFn(".repYear0Files", new.env())
  expect_identical(
    year0(rep1),
    list(cohortData = file.path(rep1, "cohortData_year0000.qs2"),
         pixelGroupMap = file.path(rep1, "pixelGroupMap_year0000.tif"))
  )

  sel <- moduleFn(".firstRepFiles", new.env())
  sam <- file.path("out", c("rep1", "rep2", "rep10"), "standAgeMap_year0701.tif")
  expect_identical(sel(sam, file.path("out", "rep1")), sam[1])
  ## a name that merely contains the folder name is not a match
  expect_identical(sel(file.path("rep1x", "a_year0701.tif"), file.path("out", "rep1")), character(0))
  ## LandWeb-style folders still work
  expect_identical(sel(file.path("o", c("rep01", "rep02"), "a.tif"), file.path("o", "rep01")),
                   file.path("o", "rep01", "a.tif"))
})
