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

## A mode = "multi" run of init and the autocorr event on stub replicates: `tables` are the replicates'
## annualSeries.csv (data.frames with a `time` column); the other files are empty stand-ins. Returns the
## messages and warnings of the run, the finished sim and the root folder.
multiRun <- function(tables, end = 600L, interval = 300L, params = list(), .env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .env)
  yrs <- unique(c(0L, moduleFn("analysesOutputsTimes", new.env())(c(0L, end), interval), end))
  for (i in seq_along(tables)) {
    d <- file.path(root, "mainSim", sprintf("rep%02d", i))
    dir.create(d, recursive = TRUE)
    utils::write.csv(tables[[i]], file.path(d, "annualSeries.csv"), row.names = FALSE)
    file.create(file.path(d, c(sprintf("cohortData_year%04d.qs2", yrs), sprintf("pixelGroupMap_year%04d.tif", yrs),
                               sprintf("standAgeMap_year%04d.tif", yrs), sprintf("vegTypeMap_year%04d.tif", yrs),
                               sprintf("flammableMap_year%04d.tif", end))))
  }
  sim <- SpaDES.core::simInit(
    times = list(start = 0, end = end),
    params = list(NRV_summary = utils::modifyList(
      list(mode = "multi", reps = seq_along(tables), summaryPeriod = c(0L, end), summaryInterval = interval,
           postprocessEvents = "autocorr", simTimes = c(0, end), timeSeriesTimes = 1:2),
      params
    )),
    modules = moduleName,
    objects = list(reportingPolygons = list(A = "placeholder")),
    paths = list(modulePath = modulePath, outputPath = file.path(root, "mainSim"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  )
  msgs <- character()
  warns <- character()
  done <- withCallingHandlers(
    SpaDES.core::spades(sim, events = list(NRV_summary = c("init", "postprocess_autocorr")), debug = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(msgs = paste(msgs, collapse = "\n"), warns = warns, sim = done, root = root)
}

## `reps` AR(1) series of `years` years from year 0, one column `x` (and `...` more columns, by name)
arTables <- function(reps, years, phi, ..., seed = 1) {
  lapply(seq_len(reps), function(r) {
    set.seed(seed * 100 + r)
    data.frame(time = seq_len(years) - 1L, x = as.numeric(stats::arima.sim(list(ar = phi), years)), ...)
  })
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
  ## development's schedule: the end-of-run events, in this order, at the last priority
  scheduleOf <- function(sim) {
    e <- as.data.frame(SpaDES.core::events(sim))[, c("eventTime", "eventType", "eventPriority")]
    e$eventTime <- as.numeric(e$eventTime)
    e$eventPriority <- as.numeric(e$eventPriority)
    rownames(e) <- NULL
    e[order(e$eventTime, e$eventPriority), ]
  }
  last <- as.numeric(SpaDES.core::.last())
  expected <- data.frame(
    eventTime = 1000,
    eventType = c("animation", "postprocess_lm", "postprocess_pm", "postprocess_lw", "postprocess", "plot"),
    eventPriority = last
  )
  expect_equal(scheduleOf(sim), expected)
  expect_identical(NROW(SpaDES.core::outputs(sim)), 0L) ## init registers nothing

  ## single mode: year 0 at the first priority + 1; the start of the summary period, the end and each
  ## timeSeriesTimes year at the last
  single <- suppressWarnings(SpaDES.core::spades(
    singleSim(end = 1000, params = list(summaryPeriod = c(700L, 1000L), summaryInterval = 50L)),
    events = list(NRV_summary = "init"), debug = FALSE
  ))
  first1 <- as.numeric(SpaDES.core::.first()) + 1
  expectedSingle <- data.frame(
    eventTime = rep(c(0, 1, 700, 1000), each = 2),
    eventType = rep(c("map_generators", "save_single"), 4),
    eventPriority = rep(c(first1, last, last, last), each = 2)
  )
  expect_equal(scheduleOf(single), expectedSingle)

  ## the files a default run registers are those it did; the annual series adds only its own
  files <- function(...) {
    done <- suppressWarnings(SpaDES.core::spades(singleSim(params = list(...)), debug = FALSE))
    basename(SpaDES.core::outputs(done)$file)
  }
  default <- files()
  expect_identical(setdiff(files(recordAnnualSeries = TRUE), default), "annualSeries.csv")
  expect_false("annualSeries.csv" %in% default)
})

test_that("no annualSeries.csv: a message, no error, and no stale csv left", {
  withr::local_package("SpaDES.core")
  sim <- multiSim(events = "autocorr")
  ## an earlier run's csv, whose replicates are gone
  stale <- file.path(dirname(SpaDES.core::outputPath(sim)), "postprocess", "csv", "autocorrelation",
                     "autocorrelation.csv")
  dir.create(dirname(stale), recursive = TRUE)
  file.create(stale)
  expect_message(
    done <- suppressWarnings(
      SpaDES.core::spades(sim, events = list(NRV_summary = c("init", "postprocess_autocorr")), debug = FALSE)
    ),
    "no annualSeries.csv found"
  )
  expect_s4_class(done, "simList")
  expect_false(file.exists(stale))
  expect_false(any(grepl("autocorrelation", SpaDES.core::outputs(done)$file)))
})

test_that("no series with a defined autocorrelation time: a message, no thinning, no -Inf", {
  const <- lapply(1:3, function(i) data.frame(time = 1:100, x = 1))
  expect_no_warning(res <- nrvtools::autocorr_summary(const, c(1, 100), 50, 10))
  expect_true(is.na(res$thin))
  expect_length(res$limiting, 1L)
  expect_true(is.na(res$limiting))

  ## the event, on a constant series
  withr::local_package("SpaDES.core")
  run <- multiRun(list(data.frame(time = 0:600, propYoung = 0.5)))
  expect_match(run$msgs, "no series has a defined autocorrelation time")
  expect_false(any(grepl("no non-missing|-Inf", run$warns)))
  expect_length(list.files(run$root, "autocorrelation", recursive = TRUE), 0L)
})

test_that("the snapshots counted by the event are the summary times of the NRV years", {
  withr::local_package("SpaDES.core")
  ## summaryPeriod 0-600 every 250, the last half: the one summary time 500 is a snapshot per rep (the count
  ## floor(300 / 250) + 1 = 2 would be wrong)
  run <- multiRun(arTables(3, 601, 0.5), interval = 250L, params = list(nrvWindow = 0.5))
  d <- utils::read.csv(list.files(run$root, "autocorrelation.csv", recursive = TRUE, full.names = TRUE))
  expect_identical(unique(d$nSnapshots), 3L)
})

test_that("the csv gives both effective snapshot counts, from nrvtools::ess_ar1(), and the snapshot count", {
  withr::local_package("SpaDES.core")
  ## summaryPeriod 0-600 every 100: 7 summary times in the NRV years of each of 4 replicates
  run <- multiRun(arTables(4, 601, 0.99), interval = 100L)
  d <- utils::read.csv(list.files(run$root, "autocorrelation.csv", recursive = TRUE, full.names = TRUE))
  nPerRep <- length(moduleFn(".nrvYears", new.env())(c(0L, 600L), 100L, 1))
  expect_identical(nPerRep, 7L)
  expect_identical(unique(d$nSnapshots), nPerRep * 4L)
  expect_equal(d$effectiveSnapshotsFinite, nrvtools::ess_ar1(nPerRep, d$r, 4))
  expect_equal(d$effectiveSnapshots, nrvtools::ess_ar1(nPerRep, d$r, 4, method = "large_sample"))
  ## nrvtools::autocorr_summary()'s columns, then the module's
  expect_identical(names(d), c("series", "phi", "tau", "ar1Misfit", "r", "rMin", "rMax", "effectiveSnapshots",
                               "effectiveSnapshotsFinite", "tauGeyer", "lagBelow0.1",
                               "label", "recommendedThin", "summaryInterval", "nSnapshots"))
  ## the message gives the finite-sample number too
  expect_match(run$msgs, paste0("about ", signif(d$effectiveSnapshots, 2), " effectively independent snapshots in the ",
                                "NRV (", signif(d$effectiveSnapshotsFinite, 2), " by the finite-sample formula)"),
               fixed = TRUE)
})

test_that("an NRV window that starts between two summary times counts the summary times inside it", {
  withr::local_package("SpaDES.core")
  ## summaryPeriod 700-1000 every 50: nrvWindow 0.55 is years 835-1000 (850, 900, 950 and 1000), 0.45 is 865-1000
  ## (900, 950 and 1000); nrvtools::autocorr_summary() stops on a snapshot outside the NRV years
  nrvYears <- moduleFn(".nrvYears", new.env())
  for (w in list(list(nrvWindow = 0.55, n = 4L), list(nrvWindow = 0.45, n = 3L))) {
    expect_length(nrvYears(c(700L, 1000L), 50L, w$nrvWindow), w$n)
    run <- multiRun(arTables(3, 1001, 0.9), end = 1000L, interval = 50L,
                    params = list(summaryPeriod = c(700L, 1000L), nrvWindow = w$nrvWindow))
    csv <- list.files(run$root, "autocorrelation.csv", recursive = TRUE, full.names = TRUE)
    expect_length(csv, 1L)
    expect_identical(unique(utils::read.csv(csv)$nSnapshots), 3L * w$n, info = w$nrvWindow)
  }
})

test_that("the event says how many replicates it read and drops a file that does not cover the NRV years", {
  withr::local_package("SpaDES.core")
  tables <- arTables(3, 601, 0.8)
  tables[[3]] <- tables[[3]][1:101, ] ## years 0-100 of 0-600
  run <- multiRun(tables)
  expect_match(run$warns, "rep03", all = FALSE)
  expect_match(run$msgs, "2 of 3 replicates")
  ## the short file did not cap the lags for the others
  d <- utils::read.csv(list.files(run$root, "autocorrelation.csv", recursive = TRUE, full.names = TRUE))
  expect_true(is.finite(d$r[1]))
})

test_that("a series missing from a replicate, or constant, is named as left out", {
  withr::local_package("SpaDES.core")
  tables <- arTables(3, 601, 0.8, flat = 2)
  tables[[1]]$lead_Pice_mar <- tables[[1]]$x
  tables[[2]]$sparse <- c(1, 2, rep(NA, 599)) ## fewer than 3 values in one replicate
  tables[[1]]$sparse <- tables[[3]]$sparse <- tables[[1]]$x
  run <- multiRun(tables)
  expect_match(run$msgs, "left out.*lead_Pice_mar")
  expect_match(run$msgs, "constant over the NRV years\\): flat")
  expect_match(run$msgs, "fewer than 3 values.*: sparse")
})

test_that("the message does not print NA when summaryInterval is longer than the NRV years", {
  withr::local_package("SpaDES.core")
  ## the NRV years of nrvWindow 0.3 (420-600) are fewer than the 300-year interval
  run <- multiRun(arTables(3, 601, 0.8), params = list(nrvWindow = 0.3))
  expect_match(run$msgs, "too short to measure the correlation")
  expect_false(grepl("r = NA|NA effectively", run$msgs))
})

test_that("the event says when the NRV is too short for the autocorrelation time it finds", {
  withr::local_package("SpaDES.core")
  ## tau about 99 over 601 years: fewer than 20 taus
  slow <- multiRun(arTables(3, 601, 0.98))
  expect_match(slow$msgs, "tend to be low")
  ## tau about 2 over 601 years
  fast <- multiRun(arTables(3, 601, 0.3))
  expect_false(grepl("tend to be low", fast$msgs))
})

test_that(".reportingCells keeps the pixels inside a reporting area smaller than the map", {
  pgm <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 400, ymin = 0, ymax = 400, vals = rep(1:4, 4))
  half <- terra::as.polygons(terra::ext(0, 200, 0, 400))
  terra::crs(half) <- terra::crs(pgm)
  cells <- moduleFn(".reportingCells", new.env())(pgm, half)
  expect_identical(sum(cells), 8L)
  expect_identical(tabulate(terra::values(pgm, mat = FALSE)[cells]), c(4L, 4L)) ## groups 1 and 2 only
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
    done <- SpaDES.core::spades(sim, events = list(NRV_summary = c("init", "postprocess_autocorr", "plot")), debug = FALSE)
  ))
  expect_match(paste(msgs, collapse = ""),
               "correlation r = .*autocorrelation time tau = .*recommended thinning [0-9]+ years")
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
