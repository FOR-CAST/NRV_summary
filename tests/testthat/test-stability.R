## The stability check asks whether a series is still trending over the last part of the NRV period.
## It is an opt-in event ('stability' in `postprocessEvents`), not in the parameter's default.
## Ported from the PE fork. Series are synthetic: 5 reps x the summary times, one metric.

## a sim with its output path in a scratch directory, and no module loaded
stabSim <- function(.env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .env)
  SpaDES.core::simInit(paths = list(
    modulePath = modulePath, outputPath = file.path(root, "mainSim"),
    cachePath = file.path(root, "cache"), inputPath = file.path(root, "in")
  ))
}

## `fn`, a module function, with P() answering `...` as the module's parameters (the module is not loaded)
withP <- function(fn, ...) {
  pars <- list(postprocessEvents = c("lm", "pm"), summaryPeriod = c(700L, 1000L), nrvWindow = 1,
               stabilityWindow = 0.5, .studyAreaName = "SA")
  pars[names(list(...))] <- list(...)
  environment(fn)$P <- function(sim, ...) pars
  fn
}

## a mode = "multi" sim of the module with its init scheduled but not run; `events` is
## postprocessEvents (NULL: the parameter's default)
initSim <- function(events = NULL, .env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .env)
  rep1 <- file.path(root, "mainSim", "rep01")
  dir.create(rep1, recursive = TRUE)
  yrs <- c(700L, 1000L)
  file.create(c(
    file.path(rep1, sprintf("cohortData_year%04d.qs2", yrs)),
    file.path(rep1, sprintf("pixelGroupMap_year%04d.tif", yrs)),
    file.path(rep1, sprintf("standAgeMap_year%04d.tif", yrs)),
    file.path(rep1, sprintf("vegTypeMap_year%04d.tif", yrs)),
    file.path(rep1, "flammableMap_year1000.tif")
  ))
  pars <- c(
    list(mode = "multi", reps = 1L, summaryPeriod = c(700L, 1000L), summaryInterval = 300L,
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

times <- seq(0, 1000, by = 100)
window <- c(500, 1000)
keys <- "metric"

makeSeries <- function(f, sd = 0, reps = 1:5, metric = "m", seed = 1) {
  set.seed(seed)
  g <- expand.grid(rep = reps, time = times)
  data.frame(metric = metric, rep = g$rep, time = g$time, value = f(g$time) + stats::rnorm(nrow(g), sd = sd))
}

trendStability <- moduleFn(".trendStability", new.env())
stabilityVerdict <- moduleFn(".stabilityVerdict", new.env())

verdict <- function(series) trendStability(series, keys, window, alpha = 0.05, minChange = 0.10)$flag

test_that("a flat series is stable", {
  expect_identical(verdict(makeSeries(function(t) 10 + 0 * t)), "stable")
})

test_that("a linear trend is still changing, with the right size", {
  s <- makeSeries(function(t) 10 + 0.01 * t, sd = 0.2)
  r <- trendStability(s, keys, window)
  expect_identical(r$flag, "still changing")
  expect_equal(r$slope, 0.01, tolerance = 0.1)
  ## the NRV range is the window's (5 to 10 here); the fitted change over 500 years is 5
  expect_equal(r$changePctRange, 100 * 5 / diff(range(s$value[s$time >= window[1]])), tolerance = 0.1)
})

test_that("a noisy flat series is stable", {
  for (seed in 1:5) {
    expect_identical(verdict(makeSeries(function(t) 10 + 0 * t, sd = 2, seed = seed)), "stable")
  }
})

test_that("the range is the window's: a transient before the window does not hide a trend in it", {
  ## values before the window run from 60 down to 10; inside it they drift from 10 to 15
  f <- function(t) 10 + ifelse(t < 500, 50 * (500 - t) / 500, 0) + 0.01 * pmax(t - 500, 0)
  s <- makeSeries(f, sd = 0.01)
  r <- trendStability(s, keys, window)
  ## the fitted change over the window (5) against the window's range (about 5), not the whole series' (about 55)
  expect_equal(r$changePctRange, 100, tolerance = 0.1)
  expect_identical(r$flag, "still changing")
})

test_that("an early transient that has levelled off does not count: only the window is tested", {
  f <- function(t) 10 + ifelse(t < 500, 50 * (500 - t) / 500, 0)
  expect_identical(verdict(makeSeries(f, sd = 0.2)), "stable")
})

test_that("fewer than 3 times in the window cannot be tested", {
  s <- makeSeries(function(t) 10 + 0.01 * t)
  r <- trendStability(s, keys, c(900, 1000))
  expect_identical(r$flag, "insufficient")
  expect_true(is.na(r$p))
})

test_that("each key combination is its own series, and NA keys are kept", {
  s <- rbind(
    makeSeries(function(t) 10 + 0 * t, metric = "flat"),
    makeSeries(function(t) 10 + 0.02 * t, metric = "rising")
  )
  s$class <- c(NA, "Pine")[1 + (s$metric == "rising")]
  r <- trendStability(s, c("metric", "class"), window)
  expect_identical(r$metric, c("flat", "rising"))
  expect_identical(r$flag, c("stable", "still changing"))
  expect_true(is.na(r$class[1]))
})

test_that("the verdict line counts the series still changing", {
  res <- data.frame(flag = c("stable", "still changing", "stable", "insufficient"))
  expect_identical(
    stabilityVerdict(res, c(850, 1000)),
    "1 of 3 series still changing; consider a longer run (1 with too few summary times to test)"
  )
  expect_identical(
    stabilityVerdict(data.frame(flag = "stable"), c(850, 1000)),
    "0 of 1 series still changing; the run is long enough for the series tested"
  )
})

test_that("the verdict does not call the run long enough when no series could be tested", {
  res <- data.frame(flag = c("insufficient", "insufficient"))
  v <- stabilityVerdict(res, c(850, 1000))
  expect_identical(
    v,
    paste0("no series could be tested: the stability window 850-1000 holds fewer than 3 summary times; ",
           "lower summaryInterval or raise stabilityWindow")
  )
  expect_false(grepl("long enough", v))
})

test_that("the stability window's summary times follow the NRV period's grid", {
  times <- moduleFn(".stabilityTimes", new.env())
  ## the module's defaults: only 900 and 1000 lie in 850-1000
  expect_identical(times(c(700L, 1000L), 100L, 1, 0.5), c(900L, 1000L))
  ## LandWeb-like settings
  expect_identical(times(c(700L, 1000L), 50L, 1, 0.5), seq(850L, 1000L, 50L))
  ## the window shrinks with nrvWindow: 925-1000
  expect_identical(times(c(700L, 1000L), 50L, 0.5, 0.5), c(950L, 1000L))
})

test_that("a stability window with fewer than 3 summary times warns, naming the window", {
  warnWindow <- moduleFn(".warnStabilityWindow", new.env())
  expect_warning(
    warnWindow(c(700L, 1000L), 100L, 1, 0.5),
    "stability window 850-1000 holds fewer than 3 summary times.*summaryInterval.*stabilityWindow"
  )
  expect_no_warning(warnWindow(c(700L, 1000L), 50L, 1, 0.5))
})

test_that("init warns, and does not stop, when stability is asked for with too short a window", {
  withr::local_package("SpaDES.core")
  run <- function(events) {
    SpaDES.core::spades(initSim(events = events), events = list(NRV_summary = "init"), debug = FALSE)
  }
  ## summaryInterval 300 over 700-1000: the window 850-1000 holds only 1000
  expect_warning(sim <- run(c("lm", "stability")), "fewer than 3 summary times")
  expect_true("postprocess_stability" %in% SpaDES.core::events(sim)$eventType)
  ## not asked for: no such warning
  w <- character(0)
  withCallingHandlers(run(c("lm", "pm")), warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("fewer than 3 summary times", w)))
})

test_that(".pages splits the series ids into pages of at most perPage", {
  pages <- moduleFn(".pages", new.env())
  p <- pages(1:210)
  expect_length(p, 18L)
  expect_identical(sort(unlist(p)), 1:210)
  expect_length(p[[18L]], 6L)
  expect_identical(pages(1:13)[[2L]], 13L)
  expect_identical(pages(character(0)), list())
})

test_that("a page draws only its own series, in a facet_wrap, with the median in bold", {
  panels <- moduleFn(".stabilityPanels", new.env())
  d <- data.frame(kind = "burn", layer = "a", poly = "a", level = "landscape", class = NA, metric = "m",
                  rep = rep(1:3, each = 2), time = rep(1:2, 3), value = c(1, 1, 1, 1, 10, 10))
  flags <- stats::setNames("stable", moduleFn(".stabilityId", new.env())(d[1L, ]))
  g <- panels(d, flags, c(1, 2), ncol = 1L)
  expect_s3_class(g$facet, "FacetWrap")
  expect_match(g$labels$caption, "bold: median")
  b <- ggplot2::ggplot_build(g)
  expect_equal(b$data[[3]]$y, c(1, 1)) ## the median of 1, 1, 10; the mean would be 4
  src <- file.path(moduleRoot, paste0(moduleName, ".R"))
  skip_if_not(file.exists(src), "module source not found")
  expect_false("facet_wrap_paginate" %in% all.names(parse(src, keep.source = FALSE)))
})

test_that("the stability series keep landscape and class rows, not one row per patch", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  m <- new.env()
  m$rptPolyNames <- "A"
  g <- expand.grid(rep = 1:2, time = c(0L, 100L))
  raw <- rbind(
    data.frame(level = "landscape", class = NA_character_, id = NA_integer_, metric = "ai", g, value = 1),
    data.frame(level = "class", class = "Pine", id = NA_integer_, metric = "ca", g, value = 2),
    data.frame(level = "patch", class = "Pine", id = 1:3, metric = "area",
               g[rep(seq_len(nrow(g)), each = 3), ], value = 3)
  )
  raw$poly <- "A"
  aggRoot <- moduleFn(".nrvAggRoot", new.env())(sim, LandWebUtils::refCodeFor("lm", "A"))
  for (r in 1:2) {
    nrvtools::write_nrv_parquet(raw[raw$rep == r, ], aggRoot, sprintf("rep%02d", r))
  }
  fn <- withP(moduleFn(".stabilitySeries", m), postprocessEvents = "lm")
  environment(fn)$open_nrv_dataset <- nrvtools::open_nrv_dataset ## attached by reqdPkgs at run time
  environment(fn)$.stabilityBurn <- function(sim) NULL ## the burn series has its own tests
  series <- fn(sim)
  expect_setequal(series$level, c("landscape", "class"))
  expect_setequal(series$kind, "lm")
  key <- c("level", "class", "metric", "rep", "time")
  expect_false(anyDuplicated(series[key]) > 0L) ## one value per rep x time per series
  expect_true(anyDuplicated(raw[key]) > 0L)
})

test_that("stability is not drawn unless postprocessEvents asks", {
  withr::local_package("SpaDES.core")
  scheduled <- function(sim) {
    sim <- suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE))
    SpaDES.core::events(sim)$eventType
  }
  expect_false("stability" %in% SpaDES.core::params(initSim())$NRV_summary$postprocessEvents)
  expect_false("postprocess_stability" %in% scheduled(initSim()))
  expect_false("postprocess_stability" %in% scheduled(initSim(events = c("lm", "pm", "lw", "am"))))
  expect_true("postprocess_stability" %in% scheduled(initSim(events = c("lm", "stability"))))
})

test_that("'stability' is an accepted postprocess event, in any case", {
  expect_no_error(moduleFn(".checkPostprocessEvents", new.env())(c("lm", "Stability")))
})

test_that("the stability figures folder is cleared only when the stability event runs", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  clear <- moduleFn(".clearFigDirs", new.env())
  f <- file.path(moduleFn(".ppFigDir", new.env())(sim, "stability", "all"), "stability_lm_p1.png")
  file.create(f)
  clear(sim, "lm")
  expect_true(file.exists(f))
  clear(sim, c("lm", "stability"))
  expect_false(file.exists(f))
})

test_that("no fireBlocks file: a message, no burn series", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  expect_message(b <- withP(moduleFn(".stabilityBurn", new.env()))(sim), "burnSummaries_fireBlocks.csv")
  expect_null(b)
})

test_that("fireBlocks beside the output path is read", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  dir.create(SpaDES.core::outputPath(sim), recursive = TRUE, showWarnings = FALSE)
  ## burnSummaries' block table has the area burned; its fire sizes have no areaBurnedHa
  blocks <- data.frame(rep = rep(1:2, each = 3), block = rep(c(2020, 2120, 2220), 2), nFires = 4L,
                       areaBurnedHa = c(100, 0, 50, 7, 7, 7))
  utils::write.csv(blocks, file.path(SpaDES.core::outputPath(sim), "burnSummaries_fireBlocks.csv"),
                   row.names = FALSE)
  burn <- withP(moduleFn(".stabilityBurn", new.env()))
  environment(burn)$.stabilityWindow <- function(sim) c(2000, 2320)
  b <- burn(sim)
  expect_setequal(b$kind, "burn")
  expect_equal(b$value[b$rep == 1], c(100, 0, 50))
  expect_equal(b$time[b$rep == 1], c(2020, 2120, 2220))
  expect_true(all(b$value[b$rep == 2] == 7))
})

test_that("a burn block that ends past the stability window is left out", {
  ## a block stamped b covers (b, b + 100]; the last one here ends at 2320
  f <- file.path(withr::local_tempdir(), "burnSummaries_fireBlocks.csv")
  blocks <- data.frame(rep = rep(1:2, each = 4), block = rep(c(1920, 2020, 2120, 2220), 2), nFires = 1L,
                       areaBurnedHa = 1)
  utils::write.csv(blocks, f, row.names = FALSE)
  series <- moduleFn(".burnSeries", new.env())
  expect_setequal(series(f, c(2000, 2320))$time, c(2020, 2120, 2220)) ## 1920 starts before the window
  expect_setequal(series(f, c(2000, 2300))$time, c(2020, 2120)) ## 2220 ends at 2320, past 2300
  expect_setequal(series(f)$time, c(1920, 2020, 2120, 2220)) ## no window: all
})

test_that("an NA among the registered files does not break the lookup", {
  blockFile <- moduleFn(".burnBlocksFile", new.env())
  other <- file.path(withr::local_tempdir(), "burnSummaries_fireBlocks.csv")
  expect_identical(blockFile(c(NA, other, NA), withr::local_tempdir()), other)
  expect_null(blockFile(c(NA_character_), withr::local_tempdir()))
})

test_that("a level filter that leaves no rows gives no series, not an error", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  m <- new.env()
  m$rptPolyNames <- "A"
  raw <- data.frame(level = "patch", class = "Pine", id = 1:3, metric = "area", rep = 1L, time = 0L,
                    value = 1, poly = "A")
  aggRoot <- moduleFn(".nrvAggRoot", new.env())(sim, LandWebUtils::refCodeFor("lm", "A"))
  nrvtools::write_nrv_parquet(raw, aggRoot, "rep01")
  fn <- withP(moduleFn(".stabilitySeries", m), postprocessEvents = "lm")
  environment(fn)$open_nrv_dataset <- nrvtools::open_nrv_dataset
  environment(fn)$.stabilityBurn <- function(sim) NULL
  expect_null(fn(sim))
})

test_that("a panel with no layer or poly name leaves them out of its label", {
  panels <- moduleFn(".stabilityPanels", new.env())
  d <- data.frame(kind = "burn", layer = NA_character_, poly = NA_character_, level = "landscape",
                  class = NA, metric = "m", rep = 1L, time = 1:2, value = 1:2)
  flags <- stats::setNames("stable", moduleFn(".stabilityId", new.env())(d[1L, ]))
  lab <- unique(panels(d, flags, c(1, 2))$data$panel)
  expect_false(grepl("NA", lab))
  expect_identical(lab, "m\nstable")
})

test_that("a registered fireBlocks file is preferred to the one beside the output path", {
  blockFile <- moduleFn(".burnBlocksFile", new.env())
  d <- withr::local_tempdir()
  beside <- file.path(d, "burnSummaries_fireBlocks.csv")
  file.create(beside)
  other <- file.path(withr::local_tempdir(), "burnSummaries_fireBlocks.csv")
  expect_identical(blockFile(c("x/burnSummaries_fireSizes_allReps.csv", other), d), other)
  expect_identical(blockFile(NULL, d), beside)
  expect_null(blockFile(character(0), withr::local_tempdir()))
})

test_that("the burn series sees a growing burn", {
  f <- file.path(withr::local_tempdir(), "burnSummaries_fireBlocks.csv")
  grow <- data.frame(rep = rep(1:5, each = 11), block = rep(seq(0, 1000, 100), 5), nFires = 1L,
                     areaBurnedHa = rep(seq(0, 1000, 100), 5) + rep(1:5, each = 11))
  utils::write.csv(grow, f, row.names = FALSE)
  g <- moduleFn(".burnSeries", new.env())(f)
  g$metric <- "areaBurnedHa per block"
  expect_identical(trendStability(g, "metric", c(500, 1000))$flag, "still changing")
})

test_that("the stability window is the last stabilityWindow of the NRV period", {
  w <- moduleFn(".stabilityWindow", new.env())
  environment(w)$.nrvPeriod <- function(sim) c(700, 1000)
  expect_equal(withP(w)(NULL), c(850, 1000))
  environment(w)$.nrvPeriod <- function(sim) c(925, 1000)
  expect_equal(withP(w, stabilityWindow = 0.2)(NULL), c(985, 1000))
})

test_that("the check writes stability.csv and one figure per page of series, each from its own series", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  ## 13 landscape series (two pages: 12 and 1) and the burn series
  g <- expand.grid(rep = 1:2, time = seq(0, 1000, 100))
  one <- function(kind, metric) {
    data.frame(kind = kind, layer = "A", poly = "A", level = "landscape", class = NA_character_,
               metric = metric, g, value = g$time / 100 + g$rep)
  }
  series <- rbind(do.call(rbind, lapply(sprintf("ai%02d", 1:13), function(m) one("lm", m))), one("burn", "areaBurnedHa per block"))
  m <- new.env()
  check <- withP(moduleFn("stabilityCheck", m), stabilityAlpha = 0.05, stabilityMinChange = 0.1)
  environment(check)$.stabilitySeries <- function(sim) series
  environment(check)$.stabilityWindow <- function(sim) c(500, 1000)
  environment(check)$.ccOverviewFigs <- function(sim) character(0) ## the overview has its own tests
  environment(check)$.metricTitle <- function(met, level) met ## the metric codes here are made up
  ## record the series each page draws
  drawn <- list()
  panels <- moduleFn(".stabilityPanels", new.env())
  idOf <- moduleFn(".stabilityId", new.env())
  environment(check)$.stabilityPanels <- function(d, ...) {
    g <- panels(d, ...)
    drawn[[length(drawn) + 1L]] <<- list(kind = unique(d$kind), ids = unique(idOf(d)), panels = unique(g$data$panel))
    g
  }
  expect_message(sim <- check(sim), "NRV_summary stability \\(500-1000\\): 14 of 14 series still changing")
  res <- utils::read.csv(file.path(moduleFn(".ppCsvDir", m)(sim, "stability", "all"), "stability.csv"))
  expect_equal(nrow(res), 14L)
  figs <- basename(list.files(moduleFn(".ppFigDir", m)(sim, "stability", "all")))
  expect_setequal(figs, c("stability_lm_p1.png", "stability_lm_p2.png", "stability_burn_p1.png", "stability_summary.png"))
  ## each page draws only its own series: 12 + 1 for lm, 1 for burn, no series twice, none missing
  expect_equal(vapply(drawn, function(p) length(p$ids), integer(1L)), c(12L, 1L, 1L))
  expect_equal(vapply(drawn, function(p) length(p$panels), integer(1L)), c(12L, 1L, 1L))
  expect_identical(vapply(drawn, function(p) p$kind, character(1L)), c("lm", "lm", "burn"))
  expect_false(anyDuplicated(unlist(lapply(drawn, `[[`, "ids"))) > 0L)
  expect_identical(sort(unlist(lapply(drawn, `[[`, "ids"))), sort(moduleFn(".stabilityId", m)(unique(series[names(series) %in% c("kind", "layer", "poly", "level", "class", "metric")]))))
  expect_true(all(file.exists(SpaDES.core::outputs(sim)$file)))
  expect_equal(length(SpaDES.core::outputs(sim)$file), 5L)
})

## the overview's synthetic envelope and current condition (also in test-current-condition-overview.R)
ccOverviewEnv <- function(polys = "ELF") {
  g <- expand.grid(poly = polys, time = c(100, 200, 300), metric = c("ai", "area_mn"), stringsAsFactors = FALSE)
  base <- ifelse(g$metric == "ai", 10, 5000)
  data.frame(poly = g$poly, time = g$time, metric = g$metric, mean = base, min = base - 2, max = base + 2)
}
ccOverviewCC <- function(polys = "ELF") {
  g <- expand.grid(poly = polys, metric = c("ai", "area_mn"), stringsAsFactors = FALSE)
  data.frame(poly = g$poly, metric = g$metric, mean = ifelse(g$metric == "ai", 11, 5000))
}

test_that("the overview pages are written under figures/stability, with the rest of the event", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  m <- new.env()
  m$rptPolyNames <- "A"
  refCode <- LandWebUtils::refCodeFor("lm", "A")
  m[[refCode]] <- ccOverviewEnv("A")
  m[[paste0(refCode, "_CC")]] <- ccOverviewCC("A")
  fn <- withP(moduleFn(".ccOverviewFigs", m), postprocessEvents = c("lm", "stability"))
  environment(fn)$.nrvPeriod <- function(sim) c(100, 300)
  files <- fn(sim)
  expect_true(length(files) >= 1L)
  expect_true(all(file.exists(files)))
  expect_true(all(grepl(file.path("figures", "stability"), files, fixed = TRUE)))
  ## no lm event, no overview
  expect_length(withP(moduleFn(".ccOverviewFigs", m), postprocessEvents = "stability")(sim), 0L)
})

test_that("stability asked for without lm says the overview needs lm", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  fn <- withP(moduleFn(".ccOverviewFigs", new.env()), postprocessEvents = "stability")
  expect_message(files <- fn(sim), "needs the lm event")
  expect_length(files, 0L)
})

## the check with its series stubbed and the real overview writer over `m`'s envelope
overviewCheck <- function(m) {
  g <- expand.grid(rep = 1:2, time = seq(0, 1000, 100))
  one <- function(metric) {
    data.frame(kind = "lm", layer = "A", poly = "A", level = "landscape", class = NA_character_,
               metric = metric, g, value = g$time / 100 + g$rep)
  }
  series <- rbind(one("ai01"), one("ai02")) ## as in the check test above
  check <- withP(moduleFn("stabilityCheck", m), postprocessEvents = c("lm", "stability"),
                 stabilityAlpha = 0.05, stabilityMinChange = 0.1)
  environment(check)$.stabilitySeries <- function(sim) series
  environment(check)$.stabilityWindow <- function(sim) c(500, 1000)
  environment(check)$.metricTitle <- function(met, level) met
  fig <- withP(moduleFn(".ccOverviewFigs", m), postprocessEvents = c("lm", "stability"))
  environment(fig)$.nrvPeriod <- function(sim) c(100, 300)
  environment(check)$.ccOverviewFigs <- fig
  check
}

test_that("the check registers the overview files as outputs", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  m <- new.env()
  m$rptPolyNames <- "A"
  refCode <- LandWebUtils::refCodeFor("lm", "A")
  m[[refCode]] <- ccOverviewEnv("A")
  m[[paste0(refCode, "_CC")]] <- ccOverviewCC("A")
  sim <- suppressMessages(overviewCheck(m)(sim))
  out <- SpaDES.core::outputs(sim)$file
  ov <- out[grepl("current_condition_A_p1.png", out, fixed = TRUE)]
  expect_length(ov, 1L)
  expect_true(file.exists(ov))
})

test_that("an overview with no comparable metric draws nothing and stability.csv is still written", {
  withr::local_package("SpaDES.core")
  sim <- stabSim()
  m <- new.env()
  m$rptPolyNames <- "A"
  refCode <- LandWebUtils::refCodeFor("lm", "A")
  m[[refCode]] <- ccOverviewEnv("A")
  cc <- ccOverviewCC("A")
  cc$metric <- "not_in_env"
  m[[paste0(refCode, "_CC")]] <- cc
  expect_no_error(sim <- suppressMessages(overviewCheck(m)(sim)))
  expect_true(file.exists(file.path(moduleFn(".ppCsvDir", m)(sim, "stability", "all"), "stability.csv")))
  expect_false(any(grepl("current_condition", SpaDES.core::outputs(sim)$file)))
})
