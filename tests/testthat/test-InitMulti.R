## InitMulti() finds each replicate's saved files in the `outputsDF` object when the run registered them
## (FireSense: folders `rep1`, `rep2`, ...) and on disk under `outputPath(sim)` when it did not
## (LandWeb: folders `rep01`, `rep02`, ..., nothing registered).

yrs <- c(700, 1000)

## the cohortData files of a replicate folder `d`, in either format: `.qs2` files (runs before 2.0.0.9046) or
## the part files of the parquet dataset grouped by year
cohortFiles <- function(d, years = yrs, format = "qs2") {
  switch(format,
         qs2 = file.path(d, sprintf("cohortData_year%04d.qs2", years)),
         parquet = file.path(d, "cohortData", sprintf("year=%04d", years), "part-0.parquet"))
}

## the files InitMulti() expects of each replicate, as paths under `d`; `cohorts` the format(s) of its cohortData
repFiles <- function(d, years = yrs, cohorts = "qs2") {
  c(
    unlist(lapply(cohorts, function(x) cohortFiles(d, years, x))),
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

## `events`: postprocessEvents (NULL: the parameter's default, lm and pm)
initMultiSim <- function(root, reps, outputsFiles = NULL, times = c(0, 1000), simTimes = c(0, 1000),
                         summaryPeriod = c(700L, 1000L), timeSeriesTimes = 701:702,
                         outputsDF = NULL, events = NULL) {
  sim <- SpaDES.core::simInit(
    times = list(start = times[1], end = times[2]),
    params = list(NRV_summary = c(
      list(mode = "multi", reps = reps, summaryPeriod = summaryPeriod, summaryInterval = 300L,
           timeSeriesTimes = timeSeriesTimes),
      if (!is.null(simTimes)) list(simTimes = simTimes),
      if (!is.null(events)) list(postprocessEvents = events)
    )),
    modules = moduleName,
    objects = c(list(reportingPolygons = list(ELF = "placeholder")),
                if (!is.null(outputsFiles)) list(outputsDF = data.frame(file = outputsFiles)),
                if (!is.null(outputsDF)) list(outputsDF = outputsDF)),
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

## init, keeping its warnings other than the one for the absent `LandTypeCC_reporting`: list(m = the module state,
## warnings =)
initWithWarnings <- function(sim) {
  w <- character(0)
  sim <- withCallingHandlers(
    SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  list(m = modOf(sim), warnings = grep("LandTypeCC_reporting", w, value = TRUE, invert = TRUE))
}

test_that("rep01 folders, no outputsDF: all reps found", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, sprintf("rep%02d", 1:3)), repFiles)))

  m <- modOf(runInit(initMultiSim(root, reps = 1:3)))
  expect_identical(m$allReps, c("1" = "rep01", "2" = "rep02", "3" = "rep03"))
  expect_identical(m$rep1Dir, c("1" = file.path(out, "rep01")))
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
  expect_identical(m$allReps, c("1" = "rep1", "2" = "rep2"))
  expect_identical(m$rep1Dir, c("1" = file.path(root, "rep1")))
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

test_that("a file listed in outputsDF but absent on disk is reported missing", {
  root <- withr::local_tempdir()
  files <- unlist(lapply(file.path(root, c("rep1", "rep2")), repFiles))
  gone <- file.path(root, "rep2", "vegTypeMap_year1000.tif")
  touchFiles(setdiff(files, gone)) ## `gone` is listed in outputsDF, but never created
  expect_error(
    runInit(initMultiSim(root, reps = 1:2, outputsFiles = files)),
    "1 simulation files appear to be missing"
  )
})

test_that("simTimes unset: year-0 and flammable rasters follow the run's start and end", {
  root <- withr::local_tempdir()
  yrs2 <- c(2720, 3020)
  files <- unlist(lapply(file.path(root, c("rep1", "rep2")), repFiles, years = yrs2))
  touchFiles(files)

  m <- modOf(runInit(initMultiSim(root, reps = 1:2, outputsFiles = files, times = c(2020, 3020),
                                  simTimes = NULL, summaryPeriod = c(2720L, 3020L),
                                  timeSeriesTimes = 2721:2722)))
  rep1 <- file.path(root, "rep1")
  expect_identical(m$fvtm0, file.path(rep1, "vegTypeMap_year2020.tif"))
  expect_identical(m$fsam0, file.path(rep1, "standAgeMap_year2020.tif"))
  expect_identical(m$ftsf0, file.path(rep1, "rstTimeSinceFire_year2020.tif"))
  expect_identical(m$flm, file.path(rep1, "flammableMap_year3020.tif"))
})

test_that("reps = 1:2 ignores rep10-rep12 on disk", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, sprintf("rep%02d", c(1:2, 10:12))), repFiles)))

  m <- modOf(runInit(initMultiSim(root, reps = 1:2)))
  expect_identical(m$allReps, c("1" = "rep01", "2" = "rep02"))
  expect_false(any(grepl("rep1[0-2]", c(m$vtm, m$sam, m$cd, m$pgm))))
  expect_length(m$vtm, 2L * length(yrs))
})

test_that("fewer reps found than asked stops", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, sprintf("rep%02d", 1:2)), repFiles)))
  expect_error(runInit(initMultiSim(root, reps = 1:3)),
               "found the files of 2 replicate\\(s\\) under .* but reps asks for 3")

  ## and with nothing at all
  root0 <- withr::local_tempdir()
  expect_error(runInit(initMultiSim(root0, reps = 1:3)),
               "found the files of 0 replicate\\(s\\) under .* but reps asks for 3")
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

test_that("aux.xml sidecars are not taken as cohort or pixelGroup files", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  files <- unlist(lapply(file.path(out, sprintf("rep%02d", 1:2)), repFiles))
  touchFiles(c(files, file.path(out, "rep01", c("pixelGroupMap_year0700.tif.aux.xml",
                                                "cohortData_year0700.qs2.aux.xml"))))

  m <- modOf(runInit(initMultiSim(root, reps = 1:2)))
  expect_false(any(grepl("xml$", c(m$cd, m$pgm))))
  expect_length(m$cd, 2L * length(yrs))
  expect_length(m$pgm, 2L * length(yrs))
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
  ## the year-0 cohorts of a parquet run, in its replicate folder
  rep1pq <- file.path(withr::local_tempdir(), "rep1")
  touchFiles(cohortFiles(rep1pq, 0, "parquet"))
  expect_identical(year0(rep1pq)$cohortData, file.path(rep1pq, "cohortData", "year=0000", "part-0.parquet"))
  ## given the cohort files InitMulti() found (here an outputsDF lists the .qs2), a stale part file beside them is not
  ## read
  rep1st <- file.path(withr::local_tempdir(), "rep1")
  touchFiles(c(cohortFiles(rep1st, 0, "qs2"), cohortFiles(rep1st, 0, "parquet")))
  expect_identical(year0(rep1st, cohortFiles(rep1st, 0, "qs2"))$cohortData, cohortFiles(rep1st, 0, "qs2"))

  sel <- moduleFn(".firstRepFiles", new.env())
  sam <- file.path("out", c("rep1", "rep2", "rep10"), "standAgeMap_year0701.tif")
  ## NA entries are dropped
  expect_identical(sel(c(NA, sam[1:2], NA), file.path("out", "rep1")), sam[1])
  ## a folder of the same name under another parent is not the first replicate's
  expect_identical(sel(file.path("other", "rep1", "a_year0701.tif"), file.path("out", "rep1")), character(0))
  expect_identical(sel(file.path("out", "sub", "rep1", "a_year0701.tif"), file.path("out", "rep1")), character(0))
  expect_identical(sel(file.path("out", "rep1", "a_year0701.tif"), file.path("out", "rep1")),
                   file.path("out", "rep1", "a_year0701.tif"))
  expect_identical(sel(sam, file.path("out", "rep1")), sam[1])
  ## a name that merely contains the folder name is not a match
  expect_identical(sel(file.path("rep1x", "a_year0701.tif"), file.path("out", "rep1")), character(0))
  ## LandWeb-style folders still work
  expect_identical(sel(file.path("o", c("rep01", "rep02"), "a.tif"), file.path("o", "rep01")),
                   file.path("o", "rep01", "a.tif"))
})

test_that("outputsDF with rows but no `file` column stops with a clear message", {
  root <- withr::local_tempdir()
  expect_error(
    runInit(initMultiSim(root, reps = 1:2, outputsDF = data.frame(path = "x/rep1/a_year0001.tif"))),
    "`outputsDF` has no `file` column",
    fixed = TRUE
  )
})

test_that(".hasOutputsDF is TRUE only for a data.frame with rows", {
  has <- moduleFn(".hasOutputsDF", new.env())
  sim <- list2env(list())
  expect_false(has(sim))
  sim$outputsDF <- data.frame(file = character(0))
  expect_false(has(sim))
  sim$outputsDF <- list(file = "a", b = 1)
  expect_false(has(sim))
  sim$outputsDF <- data.frame(file = "a")
  expect_true(has(sim))
})

test_that("files listed out of order are paired by replicate and year", {
  ## the patch metrics pair the veg-type and stand-age maps by position
  root <- withr::local_tempdir()
  files <- unlist(lapply(file.path(root, c("rep1", "rep2")), repFiles))
  touchFiles(files)
  set.seed(42)
  shuffled <- sample(files)
  expect_false(identical(shuffled, sort(files)))

  m <- modOf(runInit(initMultiSim(root, reps = 1:2, outputsFiles = shuffled)))
  byRep <- moduleFn(".filesByRep", new.env())
  yearsByRep <- function(f) lapply(byRep(f), function(x) sub(".*_year(\\d+)\\..*", "\\1", x))
  expect_identical(yearsByRep(m$vtm), list(rep1 = c("0700", "1000"), rep2 = c("0700", "1000")))
  expect_identical(yearsByRep(m$vtm), yearsByRep(m$sam))
  expect_identical(names(byRep(m$vtm)), names(byRep(m$sam)))
  expect_identical(yearsByRep(m$cd), yearsByRep(m$pgm))
})

test_that("the replicate folders are named by replicate id, and the files found are not named", {
  ## dirnamesFromSet() names its result by replicate id (PredictiveEcology/SpaDES.core#543), not by position: reps
  ## 2 and 5 are named "2" and "5", on disk and with outputsDF
  unnamed <- c("flm", "fvtm0", "fsam0", "ftsf0", "vtm", "sam", "tsf", "cd", "pgm", "allouts", "layerName")
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, c("rep02", "rep05")), repFiles)))
  m <- modOf(runInit(initMultiSim(root, reps = c(2L, 5L))))
  expect_identical(m$allReps, c("2" = "rep02", "5" = "rep05"))
  expect_identical(m$rep1Dir, c("2" = file.path(out, "rep02")))
  expect_identical(m$flm, file.path(out, "rep02", "flammableMap_year1000.tif"))
  for (x in unnamed) expect_null(names(m[[x]]), label = paste0("names(mod$", x, ")"))

  root2 <- withr::local_tempdir()
  files <- unlist(lapply(file.path(root2, c("rep2", "rep5")), repFiles))
  touchFiles(files)
  m2 <- modOf(runInit(initMultiSim(root2, reps = c(2L, 5L), outputsFiles = files)))
  expect_identical(m2$allReps, c("2" = "rep2", "5" = "rep5"))
  expect_identical(m2$rep1Dir, c("2" = file.path(root2, "rep2")))
  expect_identical(m2$flm, file.path(root2, "rep2", "flammableMap_year1000.tif"))
  for (x in unnamed) expect_null(names(m2[[x]]), label = paste0("names(mod$", x, ")"))
})

test_that("bare-number replicate folders (1/, 2/) are found through outputsDF", {
  root <- withr::local_tempdir()
  files <- unlist(lapply(file.path(root, c("1", "2")), repFiles))
  touchFiles(files)

  m <- modOf(runInit(initMultiSim(root, reps = 1:2, outputsFiles = files)))
  expect_identical(m$allReps, c("1" = "1", "2" = "2"))
  expect_identical(m$rep1Dir, c("1" = file.path(root, "1")))
  expect_length(m$vtm, 2L * length(yrs))
  expect_length(m$cd, 2L * length(yrs))
})

## ---- cohortData in either format (FOR-CAST/NRV_summary#1) ----

## the cohort file of each replicate-year, as InitMulti() selected it: "rep year format"
cdOf <- function(m) {
  part <- basename(m$cd) == "part-0.parquet"
  repDir <- ifelse(part, basename(dirname(dirname(dirname(m$cd)))), basename(dirname(m$cd)))
  year <- ifelse(part, sub("^year=", "", basename(dirname(m$cd))), sub(".*_year(\\d+)\\..*", "\\1", basename(m$cd)))
  paste(repDir, year, ifelse(part, "parquet", "qs2"))
}
pgmOf <- function(m) paste(basename(dirname(m$pgm)), sub(".*_year(\\d+)\\..*", "\\1", basename(m$pgm)))

test_that("InitMulti accepts .qs2-only replicates (old outputs)", {
  ## LandWeb's layout: rep01 folders, nothing registered
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, sprintf("rep%02d", 1:3)), repFiles)))

  r <- initWithWarnings(initMultiSim(root, reps = 1:3))
  m <- r$m
  expect_false(any(grepl("both formats", r$warnings)))
  expect_identical(m$cd, unlist(lapply(file.path(out, sprintf("rep%02d", 1:3)), cohortFiles)))
  expect_identical(m$pgm, file.path(out, rep(sprintf("rep%02d", 1:3), each = 2L),
                                    sprintf("pixelGroupMap_year%04d.tif", yrs)))
})

test_that("InitMulti finds parquet cohorts on disk and in outputsDF, one per replicate and year", {
  ## rep01 saved parquet, rep02 .qs2, rep03 both: rep03's parquet files are taken, and each cohort file is in the
  ## place of the pixelGroupMap of its replicate and year (BC seral pairs them by position)
  layout <- list(rep01 = "parquet", rep02 = "qs2", rep03 = c("qs2", "parquet"))
  want <- c(paste(rep(names(layout), each = 2L), c("0700", "1000"), rep(c("parquet", "qs2", "parquet"), each = 2L)))
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  written <- unlist(Map(function(r, f) repFiles(file.path(out, r), cohorts = f), names(layout), layout))
  touchFiles(written)
  r <- initWithWarnings(initMultiSim(root, reps = 1:3))
  m <- r$m
  expect_identical(cdOf(m), want)
  expect_identical(sub(" [a-z0-9]+$", "", cdOf(m)), pgmOf(m))
  expect_setequal(m$cohortFiles, grep("cohortData", written, value = TRUE))
  ## rep03 has both formats, at two years: one warning, which names the parquet files read
  both <- grep("both formats", r$warnings, value = TRUE)
  expect_length(both, 1L)
  expect_match(both, "2 replicate-year(s) with cohortData in both formats", fixed = TRUE)
  expect_match(both, file.path(out, "rep03", "cohortData", "year=0700", "part-0.parquet"), fixed = TRUE)

  ## the same from outputsDF, listed out of order
  layout2 <- stats::setNames(layout, c("rep1", "rep2", "rep3"))
  root2 <- withr::local_tempdir()
  files <- unlist(Map(function(r, f) repFiles(file.path(root2, r), cohorts = f), names(layout2), layout2))
  touchFiles(files)
  set.seed(1)
  r2 <- initWithWarnings(initMultiSim(root2, reps = 1:3, outputsFiles = sample(files)))
  m2 <- r2$m
  expect_identical(cdOf(m2), sub("^rep0", "rep", want))
  expect_identical(sub(" [a-z0-9]+$", "", cdOf(m2)), pgmOf(m2))
  expect_identical(m2$allReps, c("1" = "rep1", "2" = "rep2", "3" = "rep3"))
  expect_identical(m2$rep1Dir, c("1" = file.path(root2, "rep1")))
  expect_match(grep("both formats", r2$warnings, value = TRUE), "2 replicate-year(s)", fixed = TRUE)
})

test_that("a replicate-year with cohorts in neither format is missing, and the message names both", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  files <- unlist(lapply(file.path(out, sprintf("rep%02d", 1:2)), repFiles, cohorts = "parquet"))
  gone <- cohortFiles(file.path(out, "rep02"), 1000, "parquet")
  touchFiles(setdiff(files, gone))
  rep2 <- file.path(out, "rep02")
  expect_error(
    runInit(initMultiSim(root, reps = 1:2)),
    paste0("1 simulation files appear to be missing:\n", file.path(rep2, "cohortData_year1000.qs2"), " or ",
           file.path(rep2, "cohortData", "year=1000", "part-0.parquet"), "$")
  )
  ## nor is a part file found when the outputsDF does not list it
  root2 <- withr::local_tempdir()
  files2 <- unlist(lapply(file.path(root2, c("rep1", "rep2")), repFiles, cohorts = "parquet"))
  touchFiles(files2)
  listed <- setdiff(files2, cohortFiles(file.path(root2, "rep2"), 1000, "parquet"))
  expect_error(
    runInit(initMultiSim(root2, reps = 1:2, outputsFiles = listed)),
    "1 simulation files appear to be missing:.*rep2/cohortData_year1000.qs2 or .*rep2/cohortData/year=1000/part-0"
  )
})

test_that(".discoverOutputFiles adds the cohort part files to a LandWeb-like listing, and nothing else", {
  discover <- moduleFn(".discoverOutputFiles", new.env())
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  ## what a LandWeb mainSim folder holds: replicate folders with figure and dataset subfolders, files at the top
  reps <- file.path(out, c("rep01", "rep02"))
  touchFiles(c(
    unlist(lapply(reps, repFiles)),
    file.path(out, c("burnSummaries_fireSizes_allReps.csv", "meanAnnualCumulBurnMap.tif")),
    file.path(out, "_cc", "rep01", "vegTypeMap_year0000.tif"),
    file.path(out, "figures", "NRV_summary", "a.png"),
    file.path(reps[1], "burnSummaries_fireSizes", "replicate=1", "part-0.parquet"),
    file.path(reps[1], "vegetation-transitions", "zone=A", "time=700", "part-0.parquet"),
    ## named like cohort part files, but not in a replicate's own cohortData folder
    file.path(reps[1], "figures", "cohortData", "year=0700", "part-0.parquet"),
    file.path(out, "cohortData", "year=0700", "part-0.parquet")
  ))
  sim <- initMultiSim(root, reps = 1:2)
  before <- as.character(fs::dir_ls(out, recurse = 1, type = "file")) ## the listing of 2.0.0.9045
  expect_identical(discover(sim), before)

  parts <- c(cohortFiles(reps[1], yrs, "parquet"), cohortFiles(reps[2], 700, "parquet"))
  touchFiles(c(parts, file.path(reps[1], "cohortData", c("year=0700/part-1.parquet", "year=07x0/part-0.parquet",
                                                        "part-0.parquet"))))
  expect_identical(discover(sim), c(before, parts))
})

test_that("reps = 1:2 ignores rep10-rep12 parquet cohorts on disk", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  touchFiles(unlist(lapply(file.path(out, sprintf("rep%02d", c(1:2, 10:12))), repFiles, cohorts = "parquet")))

  m <- modOf(runInit(initMultiSim(root, reps = 1:2)))
  expect_identical(m$allReps, c("1" = "rep01", "2" = "rep02"))
  rep12 <- file.path(out, c("rep01", "rep02"))
  expect_identical(m$cd, unlist(lapply(rep12, cohortFiles, format = "parquet")))
  expect_setequal(m$cohortFiles, unlist(lapply(rep12, cohortFiles, format = "parquet")))
})

test_that("a pixelGroupMap with no cohort file of its folder and year stops InitMulti, naming it, if bc or st runs", {
  root <- withr::local_tempdir()
  files <- unlist(lapply(file.path(root, c("rep1", "rep2")), repFiles, cohorts = "parquet"))
  extra <- file.path(root, "rep1", "sub", "pixelGroupMap_year0700.tif")
  touchFiles(c(files, extra))
  ## only the bc and st events read the cohort files and pixel group maps: lm and pm (the default) run as before
  expect_no_error(runInit(initMultiSim(root, reps = 1:2, outputsFiles = c(files, extra))))
  expect_no_error(runInit(initMultiSim(root, reps = 1:2, outputsFiles = c(files, extra), events = c("lm", "lw", "am"))))
  for (ev in c("st", "bc", "BC")) {
    expect_error(
      runInit(initMultiSim(root, reps = 1:2, outputsFiles = c(files, extra), events = c("lm", ev))),
      paste0("no cohortData file in the folder and of the year of these pixelGroupMap files:\n", extra, "$"),
      info = ev
    )
  }
})

test_that("a folder named like the pixelGroupMap files does not make other files pixel group maps", {
  root <- withr::local_tempdir(pattern = "pixelGroupMap_")
  out <- file.path(root, "mainSim")
  reps <- file.path(out, c("rep01", "rep02"))
  touchFiles(unlist(lapply(reps, repFiles)))

  m <- modOf(runInit(initMultiSim(root, reps = 1:2)))
  expect_identical(m$pgm, file.path(rep(reps, each = 2L), sprintf("pixelGroupMap_year%04d.tif", yrs)))
  expect_identical(m$cd, unlist(lapply(reps, cohortFiles)))
})

test_that("a part file the outputsDF lists but that is not on disk is not taken: the .qs2 beside it is", {
  root <- withr::local_tempdir()
  files <- unlist(lapply(file.path(root, c("rep1", "rep2")), repFiles, cohorts = c("qs2", "parquet")))
  touchFiles(files)
  gone <- cohortFiles(file.path(root, "rep1"), 700, "parquet")
  file.remove(gone)
  r <- initWithWarnings(initMultiSim(root, reps = 1:2, outputsFiles = files, events = c("lm", "st")))
  expect_identical(cdOf(r$m), c("rep1 0700 qs2", "rep1 1000 parquet", "rep2 0700 parquet", "rep2 1000 parquet"))
  expect_false(gone %in% r$m$cohortFiles)
  ## the other three replicate-years have both formats
  expect_match(grep("both formats", r$warnings, value = TRUE), "3 replicate-year(s)", fixed = TRUE)
})

test_that(".cohortFile pairs the files of a folder and year, whatever the year's padding", {
  cohortFile <- moduleFn(".cohortFile", new.env())
  rep1 <- file.path("o", "rep01")
  rep2 <- file.path("o", "rep02")
  built <- "cohortData_year0700.qs2"
  padded <- file.path(rep1, "cohortData_year700.qs2")
  ## a file padded differently beside the pixelGroupMap of that year
  expect_identical(cohortFile(rep1, 700L, built, files = padded), padded)
  ## its replicate folder counts: another folder's file of that year is not taken
  other <- file.path(rep2, "cohortData_year700.qs2")
  expect_identical(cohortFile(rep1, 700L, built, files = other), file.path(rep1, built))
  expect_identical(cohortFile(c(rep1, rep2), c(700L, 700L), c(built, built), files = c(padded, other)),
                   c(padded, other))
  ## another year's is not either
  expect_identical(cohortFile(rep1, 1000L, "cohortData_year1000.qs2", files = padded),
                   file.path(rep1, "cohortData_year1000.qs2"))
  ## the file named as the caller names it is kept when it is among them
  ## (the other candidate sorts before it, so taking the first by name would give that one)
  exact <- file.path(rep1, built)
  early <- file.path(rep1, "cohortData_year00700.qs2")
  expect_true(early < exact)
  expect_identical(cohortFile(rep1, 700L, built, files = c(early, exact)), exact)
  ## none named so: the first by name
  expect_identical(cohortFile(rep1, 700L, built, files = c(padded, early)), early)
  ## without candidates, the path is built as before
  expect_identical(cohortFile(rep1, 700L, built), file.path(rep1, built))
  ## the part file wins over a .qs2 of any padding, when it exists and is a candidate
  root <- withr::local_tempdir()
  r1 <- file.path(root, "rep01")
  pq <- touchFiles(cohortFiles(r1, 700, "parquet"))
  qs <- touchFiles(file.path(r1, "cohortData_year700.qs2"))
  expect_identical(cohortFile(r1, 700L, built, files = c(qs, pq)), pq)
  expect_identical(cohortFile(r1, 700L, built, files = qs), qs)
})

test_that("InitMulti pairs a cohort file padded differently from its pixelGroupMap", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  reps <- file.path(out, sprintf("rep%02d", 1:2))
  files <- unlist(lapply(reps, repFiles))
  files <- sub("cohortData_year0700.qs2", "cohortData_year700.qs2", files, fixed = TRUE)
  touchFiles(files)
  m <- modOf(runInit(initMultiSim(root, reps = 1:2, events = c("lm", "pm", "bc"))))
  expect_identical(m$cd, file.path(rep(reps, each = 2L), c("cohortData_year700.qs2", "cohortData_year1000.qs2")))
  fileYear <- moduleFn(".fileYear", new.env())
  expect_identical(fileYear(m$cd), fileYear(m$pgm))
})

test_that("InitMulti warns of several .qs2 cohort files of one replicate-year, and names the one used", {
  root <- withr::local_tempdir()
  out <- file.path(root, "mainSim")
  rep1 <- file.path(out, "rep01")
  touchFiles(c(repFiles(rep1), file.path(rep1, "cohortData_year700.qs2")))
  r <- initWithWarnings(initMultiSim(root, reps = 1L))
  w <- grep("several .qs2", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "cohortData_year700.qs2", fixed = TRUE)
  expect_match(w, paste0("used:\n", file.path(rep1, "cohortData_year0700.qs2")), fixed = TRUE)
  expect_identical(cdOf(r$m)[1], "rep01 0700 qs2")
  ## one file per replicate-year: no warning
  root2 <- withr::local_tempdir()
  touchFiles(repFiles(file.path(root2, "mainSim", "rep01")))
  expect_false(any(grepl("several .qs2", initWithWarnings(initMultiSim(root2, reps = 1L))$warnings)))
})
