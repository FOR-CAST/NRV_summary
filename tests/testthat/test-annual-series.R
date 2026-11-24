## The annual landscape series (`recordAnnualSeries`) and the autocorrelation check (the `autocorr`
## event). Both are opt-in: with the defaults a run schedules neither. Ported from the PE fork (PE#14).

## a mode = "multi" sim of the module with its init scheduled but not run. The replicates' files are
## empty stand-ins for each year in `yrs`; `events` is postprocessEvents (NULL: the parameter's default).
multiSim <- function(events = NULL, yrs = c(700L, 1000L), interval = 300L, period = c(700L, 1000L),
                     .env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .env)
  rep1 <- file.path(root, "mainSim", "rep01")
  dir.create(rep1, recursive = TRUE)
  file.create(c(
    file.path(rep1, sprintf("cohortData_year%04d.qs2", yrs)),
    file.path(rep1, sprintf("pixelGroupMap_year%04d.tif", yrs)),
    file.path(rep1, sprintf("standAgeMap_year%04d.tif", yrs)),
    file.path(rep1, sprintf("vegTypeMap_year%04d.tif", yrs)),
    file.path(rep1, "flammableMap_year1000.tif")
  ))
  pars <- c(
    list(mode = "multi", reps = 1L, summaryPeriod = period, summaryInterval = interval,
         timeSeriesTimes = 701:702, simTimes = c(0, 1000)),
    if (!is.null(events)) list(postprocessEvents = events)
  )
  SpaDES.core::simInit(
    times = list(start = 0, end = 1000),
    params = list(NRV_summary = pars),
    modules = moduleName,
    objects = list(reportingPolygons = list(A = "placeholder")),
    paths = list(modulePath = modulePath, outputPath = file.path(root, "mainSim"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  )
}

## a mode = "single" sim of the small fixture landscape, 0 to `end`
singleSim <- function(end = 4, params = list(), out = withr::local_tempdir(.local_envir = .env),
                      .env = parent.frame()) {
  SpaDES.core::simInit(
    times = list(start = 0, end = end),
    params = list(NRV_summary = utils::modifyList(
      list(mode = "single", summaryPeriod = c(0L, as.integer(end)), summaryInterval = 2L,
           timeSeriesTimes = 1, ccStandAgeSource = "cohortData"),
      params
    )),
    modules = moduleName,
    objects = singleRunFixture()$objects,
    paths = utils::modifyList(testPaths, list(outputPath = out))
  )
}

initEvents <- function(sim) {
  sim <- suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE))
  SpaDES.core::events(sim)$eventType
}

test_that("recordAnnualSeries defaults to FALSE: no annual_series events", {
  withr::local_package("SpaDES.core")
  expect_false(SpaDES.core::params(singleSim())$NRV_summary$recordAnnualSeries)
  expect_false("annual_series" %in% initEvents(singleSim()))
  expect_true("annual_series" %in% initEvents(singleSim(params = list(recordAnnualSeries = TRUE))))

  ## and a whole default run writes no annualSeries.csv
  out <- withr::local_tempdir()
  sim <- suppressWarnings(SpaDES.core::spades(singleSim(out = out), debug = FALSE))
  expect_false(file.exists(file.path(out, "annualSeries.csv")))
  expect_false("annual_series" %in% as.data.frame(SpaDES.core::completed(sim))$eventType)
  expect_false(any(basename(SpaDES.core::outputs(sim)$file) == "annualSeries.csv"))
})

test_that("single mode writes one annualSeries.csv row per year and registers it", {
  out <- withr::local_tempdir()
  sim <- suppressWarnings(
    SpaDES.core::spades(singleSim(out = out, params = list(recordAnnualSeries = TRUE)), debug = FALSE)
  )
  f <- file.path(out, "annualSeries.csv")
  expect_true(file.exists(f))
  d <- utils::read.csv(f)
  expect_identical(d$time, 0:4)
  expect_true(all(c("propYoung", "propOld", "meanStandAge", "totalBiomassTg", "lead_Pice_gla") %in% names(d)))
  expect_true(any(basename(SpaDES.core::outputs(sim)$file) == "annualSeries.csv"))
})

test_that("the autocorr event is not run unless postprocessEvents asks", {
  withr::local_package("SpaDES.core")
  expect_false("autocorr" %in% SpaDES.core::params(multiSim())$NRV_summary$postprocessEvents)
  expect_false("postprocess_autocorr" %in% initEvents(multiSim()))
  expect_false("postprocess_autocorr" %in% initEvents(multiSim(events = c("lm", "pm", "lw", "am"))))
  expect_true("postprocess_autocorr" %in% initEvents(multiSim(events = c("lm", "autocorr"))))
  expect_true("postprocess_autocorr" %in% initEvents(multiSim(events = "Autocorr")))
})

test_that("'autocorr' is an accepted postprocess event, in any case", {
  expect_no_error(moduleFn(".checkPostprocessEvents", new.env())(c("lm", "Autocorr")))
})

test_that("with LandWeb's parameters, init schedules the events it did before and registers no new files", {
  withr::local_package("SpaDES.core")
  yrs <- seq(700L, 1000L, 50L)
  events <- c("lm", "pm", "lw", "am")
  sim <- suppressWarnings(SpaDES.core::spades(
    multiSim(events = events, yrs = yrs, interval = 50L), events = list(NRV_summary = "init"), debug = FALSE
  ))
  expect_setequal(
    SpaDES.core::events(sim)$eventType,
    c("animation", "postprocess_lm", "postprocess_pm", "postprocess_lw", "postprocess", "plot")
  )
  expect_false(any(grepl("annualSeries|autocorrelation", SpaDES.core::outputs(sim)$file)))

  single <- singleSim(end = 10, params = list(summaryPeriod = c(5L, 10L), summaryInterval = 5L))
  expect_setequal(initEvents(single), c("map_generators", "save_single"))
})

test_that("no annualSeries.csv: a message, no error", {
  withr::local_package("SpaDES.core")
  sim <- suppressWarnings(SpaDES.core::spades(
    multiSim(events = "autocorr"), events = list(NRV_summary = "init"), debug = FALSE
  ))
  modEnv <- new.env()
  modEnv$annualSeriesFiles <- character(0)
  autocorrSummary <- moduleFn("autocorrSummary", modEnv)
  expect_message(res <- autocorrSummary(sim), "no annualSeries.csv found")
  expect_s4_class(res, "simList")
})

test_that("multi mode reads the reps' annual series and writes the autocorrelation csv and figure", {
  withr::local_package("SpaDES.core")
  root <- withr::local_tempdir()
  end <- 600L
  period <- c(0L, end)
  yrs <- c(0L, 300L, 600L)
  files <- unlist(lapply(1:3, function(r) {
    d <- file.path(root, "mainSim", sprintf("rep%02d", r))
    dir.create(d, recursive = TRUE)
    set.seed(r)
    ts <- data.frame(time = 0:end, propYoung = as.numeric(stats::arima.sim(list(ar = 0.9), end + 1)),
                     lead_Pice_gla = as.numeric(stats::arima.sim(list(ar = 0.3), end + 1)))
    utils::write.csv(ts, file.path(d, "annualSeries.csv"), row.names = FALSE)
    f <- file.path(d, c(sprintf("cohortData_year%04d.qs2", yrs), sprintf("pixelGroupMap_year%04d.tif", yrs),
                        sprintf("standAgeMap_year%04d.tif", yrs), sprintf("vegTypeMap_year%04d.tif", yrs),
                        "flammableMap_year0600.tif"))
    file.create(f)
    c(file.path(d, "annualSeries.csv"), f)
  }))
  sim <- SpaDES.core::simInit(
    times = list(start = 0, end = end),
    params = list(NRV_summary = list(mode = "multi", reps = 1:3, summaryPeriod = period, summaryInterval = 300L,
                                     nrvWindow = 1, postprocessEvents = "autocorr", simTimes = c(0, end),
                                     timeSeriesTimes = 1:2, .studyAreaName = "stub")),
    modules = moduleName,
    objects = list(outputsDF = data.table::data.table(file = files, saveTime = end),
                   reportingPolygons = list(ELF = "x"), sppEquiv = singleRunFixture()$objects$sppEquiv),
    paths = list(modulePath = modulePath, outputPath = file.path(root, "mainSim"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  )
  ## a figure of an earlier run is cleared at the start of the event
  stale <- file.path(root, "postprocess", "figures", "autocorrelation", "all", "old.png")
  dir.create(dirname(stale), recursive = TRUE)
  file.create(stale)

  msgs <- testthat::capture_messages(suppressWarnings(
    done <- SpaDES.core::spades(sim, events = list(NRV_summary = c("init", "postprocess_autocorr")), debug = FALSE)
  ))
  expect_match(paste(msgs, collapse = ""),
               "correlation r = .*fitted autocorrelation time tau = .*recommended thinning [0-9]+ years")
  csv <- list.files(root, "autocorrelation.csv", recursive = TRUE, full.names = TRUE)
  expect_length(csv, 1L)
  d <- utils::read.csv(csv)
  expect_true(all(c("series", "phi", "tau", "tauGeyer", "ar1Misfit", "r", "rMin", "rMax", "label",
                    "recommendedThin") %in% names(d)))
  expect_true("Proportion of forest young" %in% d$label)
  expect_gt(d$recommendedThin[1], 5)
  png <- list.files(root, "autocorrelation.png", recursive = TRUE, full.names = TRUE)
  expect_length(png, 1L)
  expect_false(file.exists(stale))

  ## both files are registered as outputs
  registered <- SpaDES.core::outputs(done)$file
  expect_true(all(c("autocorrelation.csv", "autocorrelation.png") %in% basename(registered)))
})

test_that("the autocorrelation figures folder is cleared only when the autocorr event is among those run", {
  withr::local_package("SpaDES.core")
  sim <- SpaDES.core::simInit(paths = testPaths)
  clear <- moduleFn(".clearFigDirs", new.env())
  f <- file.path(moduleFn(".ppFigDir", new.env())(sim, "autocorrelation", "all"), "autocorrelation.png")
  file.create(f)
  clear(sim, c("lm", "stability"))
  expect_true(file.exists(f))
  clear(sim, "autocorr")
  expect_false(file.exists(f))
})

test_that("the plot event's clearing leaves the autocorrelation figures, which the earlier event drew", {
  ## plotFun() clears before the plot event draws, after the autocorr event wrote its figure
  src <- file.path(moduleRoot, paste0(moduleName, ".R"))
  skip_if_not(file.exists(src), "module source not found")
  code <- paste(deparse(moduleFn("plotFun", new.env())), collapse = "\n")
  expect_match(code, 'setdiff\\(events, "autocorr"\\)')
})
