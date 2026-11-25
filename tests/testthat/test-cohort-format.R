## Mode = "single" saves cohortData at each summary year as one `cohortData_year<YYYY>.qs2` file per year, by default
## and as every run before 2.0.0.9046 did, or with `cohortDataFormat = "parquet"` as a parquet dataset grouped by
## simulation time, `cohortData/year=<YYYY>/part-0.parquet` (FOR-CAST/NRV_summary#1). Mode = "multi" reads both (see
## test-InitMulti.R, test-reporting-groups.R and test-standMetrics.R).

## a mode = "single" run of the shared fixture from 0 to 1000, saving at 0, 700 and 1000 (and 701 for the time
## series, which saves no cohorts); `format` is cohortDataFormat (NULL: the parameter's default)
runSave <- function(format = NULL, .env = parent.frame()) {
  out <- withr::local_tempdir(.local_envir = .env)
  fx <- singleRunFixture()
  sim <- SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 1000),
    params = list(NRV_summary = c(
      list(mode = "single", summaryPeriod = c(700L, 1000L), summaryInterval = 300L, timeSeriesTimes = 701,
           ccStandAgeSource = "cohortData"),
      if (!is.null(format)) list(cohortDataFormat = format)
    )),
    modules = moduleName,
    objects = fx$objects,
    paths = utils::modifyList(testPaths, list(outputPath = out))
  )
  list(out = out, cohortData = fx$cohortData,
       registered = normalizePath(SpaDES.core::outputs(sim)$file, winslash = "/", mustWork = FALSE))
}

test_that("save_single writes cohortData/year=0700/part-0.parquet and registers it", {
  r <- runSave("parquet")
  f <- file.path(r$out, "cohortData", "year=0700", "part-0.parquet")
  expect_true(file.exists(f))
  expect_true(normalizePath(f, winslash = "/") %in% r$registered)
  expect_equal(as.data.frame(nrvtools::read_cohort_data(f)), as.data.frame(r$cohortData))
  ## the cohorts of every summary year, year 0 among them, and no .qs2 cohorts
  expect_identical(list.files(file.path(r$out, "cohortData")), c("year=0000", "year=0700", "year=1000"))
  expect_identical(sum(basename(r$registered) == "part-0.parquet"), 3L)
  expect_length(list.files(r$out, "^cohortData_year"), 0L)
})

test_that("by default save_single writes and registers cohortData_year0700.qs2 as before", {
  qs <- runSave() ## cohortDataFormat = "qs2", the default
  f <- file.path(qs$out, "cohortData_year0700.qs2")
  expect_true(file.exists(f))
  expect_true(normalizePath(f, winslash = "/") %in% qs$registered)
  expect_identical(qs2::qs_read(f)$B, qs$cohortData$B)
  expect_false(dir.exists(file.path(qs$out, "cohortData")))
  ## the files of 2.0.0.9045, which "qs2" given explicitly writes too; "parquet" changes the cohort files only
  expect_setequal(basename(qs$registered), c(
    "sppColorVect_year0000.qs2", "sppEquiv_year0000.qs2", "speciesLayers_year0000.tif", "flammableMap_year1000.tif",
    sprintf("cohortData_year%04d.qs2", c(0, 700, 1000)), sprintf("pixelGroupMap_year%04d.tif", c(0, 700, 1000)),
    sprintf("standAgeMap_year%04d.tif", c(0, 700, 701, 1000)), sprintf("vegTypeMap_year%04d.tif", c(0, 700, 701, 1000))
  ))
  rel <- function(r) substring(r$registered, nchar(normalizePath(r$out, winslash = "/")) + 2L)
  expect_identical(rel(runSave("qs2")), rel(qs))
  pq <- runSave("parquet")
  expect_identical(setdiff(rel(qs), rel(pq)), sprintf("cohortData_year%04d.qs2", c(0, 700, 1000)))
  expect_identical(setdiff(rel(pq), rel(qs)), sprintf("cohortData/year=%04d/part-0.parquet", c(0, 700, 1000)))
})

test_that("init stops on a cohortDataFormat other than parquet or qs2, in either mode", {
  withr::local_package("SpaDES.core")
  initWith <- function(mode) {
    sim <- SpaDES.core::simInit(
      times = list(start = 0, end = 1000),
      params = list(NRV_summary = list(mode = mode, cohortDataFormat = "rds", summaryPeriod = c(700L, 1000L),
                                       timeSeriesTimes = 701)),
      modules = moduleName,
      objects = list(reportingPolygons = list(A = "placeholder")),
      paths = utils::modifyList(testPaths, list(outputPath = withr::local_tempdir()))
    )
    suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE))
  }
  for (mode in c("single", "multi")) {
    expect_error(initWith(mode), "cohortDataFormat must be \"parquet\" or \"qs2\"", fixed = TRUE, info = mode)
  }
})
