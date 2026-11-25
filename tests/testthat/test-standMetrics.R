## The stand-level measures of the opt-in `st` event, from synthetic snapshot values, and the NRV comparison
## tables that go with it. Ported from the PE fork (PE#19). Both are opt-in: with LandWeb's parameters a run
## schedules no `st` event and its lm event writes no comparison table.

standRow <- function(...) moduleFn(".standRow", new.env())(...)
standLabel <- function(...) moduleFn(".standLabel", new.env())(...)
standClass <- function(...) moduleFn(".standClass", new.env())(...)
standSeriesPanels <- function(...) moduleFn(".standSeriesPanels", new.env())(...)
nrvComparison <- function(...) moduleFn(".nrvComparison", new.env())(...)
positionColour <- function(...) moduleFn(".positionColour", new.env())(...)

## a sim with its output path in a scratch directory, and no module loaded
stSim <- function(.env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .env)
  SpaDES.core::simInit(paths = list(
    modulePath = modulePath, outputPath = file.path(root, "mainSim"),
    cachePath = file.path(root, "cache"), inputPath = file.path(root, "in")
  ))
}

## `fn`, a module function, with P() answering `...` as the module's parameters (the module is not loaded)
stP <- function(fn, ...) {
  pars <- list(postprocessEvents = c("lm", "pm"), summaryPeriod = c(700L, 1000L), nrvWindow = 1,
               .studyAreaName = "SA")
  pars[names(list(...))] <- list(...)
  environment(fn)$P <- function(sim, ...) pars
  fn
}

test_that("a snapshot's stand-level row: age classes of forest, forested and leading shares of the area, biomass", {
  ## 6 reporting-area pixels: 4 forested (ages 10, 50, 90, 150), one with no age, one with pixel group 0
  sam <- c(10, 50, 90, 150, NA, 30)
  pgm <- c(1L, 1L, 2L, 2L, 3L, 0L)
  vtm <- c("Pice_mar", "Pice_mar", "Mixed", "Pice_mar", NA, NA)
  cd <- data.frame(pixelGroup = c(1L, 1L, 2L), B = c(100, 200, 1000)) ## g/m2: 300 in group 1, 1000 in group 2
  cls <- c("Young", "Immature", "Mature", "Old")
  r <- standRow(sam, vtm, pgm, cd, pixelHa = 6.25, ageClassCutOffs = c(0, 40, 80, 120), ageClasses = cls)
  v <- setNames(r$value, r$metric)
  expect_equal(v[paste0("age_", cls)], setNames(rep(25, 4), paste0("age_", cls)))
  expect_equal(v[["forested"]], 100 * 4 / 6)
  expect_equal(v[["lead_Pice_mar"]], 100 * 3 / 6)
  expect_equal(v[["lead_Mixed"]], 100 * 1 / 6)
  ## B over the pixels: 2 x 300 + 2 x 1000 g/m2, x 6.25 ha x 1e4 m2/ha, in Tg
  expect_equal(v[["biomassTg"]], (2 * 300 + 2 * 1000) * 6.25 * 1e4 / 1e12)
})

test_that("current biomass excludes a cell masked as urban", {
  ## the current condition's stand age is masked where urban (NA), its pixel group map is not: the cell has
  ## cohorts (pixel group 3, 5000 g/m2) but is not forest, so its biomass is not counted
  sam <- c(10, 50, NA)
  pgm <- c(1L, 1L, 3L)
  vtm <- c("Pice_mar", "Pice_mar", NA)
  cd <- data.frame(pixelGroup = c(1L, 3L), B = c(300, 5000))
  r <- standRow(sam, vtm, pgm, cd, pixelHa = 6.25, ageClassCutOffs = c(0, 40), ageClasses = c("Young", "Old"))
  v <- setNames(r$value, r$metric)
  expect_equal(v[["biomassTg"]], 2 * 300 * 6.25 * 1e4 / 1e12)
  expect_equal(v[["forested"]], 100 * 2 / 3)
})

test_that("snapshot files pair by rep folder and year", {
  f <- function(kind, ext) c(sprintf("o/rep2/%s_year2120.%s", kind, ext), sprintf("o/rep1/%s_year2020.%s", kind, ext),
                             sprintf("o/rep1/%s_year2120.%s", kind, ext))
  s <- moduleFn(".snapshotTable", new.env())(f("standAgeMap", "tif"), f("vegTypeMap", "tif"), f("pixelGroupMap", "tif"),
                                             f("cohortData", "qs2"))
  expect_identical(nrow(s), 3L)
  i <- s$rep == "rep2"
  expect_identical(s$time[i], 2120)
  expect_identical(s$cd[i], "o/rep2/cohortData_year2120.qs2")
  expect_identical(s$vtm[i], "o/rep2/vegTypeMap_year2120.tif")
})

test_that("a type absent from a snapshot leads 0% of it, also one found only in the current condition", {
  rows <- data.frame(rep = c(1, 1, 2), time = c(0, 100, 0), metric = c("lead_A", "lead_B", "lead_A"), value = c(5, 7, 3))
  out <- moduleFn(".fillAbsentMetrics", new.env())(rows, "lead_C")
  expect_identical(nrow(out), 2L * 2L * 3L)
  expect_equal(out$value[out$rep == 2 & out$time == 100 & out$metric == "lead_A"], 0)
  expect_equal(sum(out$value[out$metric == "lead_C"]), 0)
})

test_that("the current condition is inside, above or below the NRV range of the NRV years only", {
  env <- data.frame(poly = "E", metric = rep(c("a", "b", "c"), each = 3), time = rep(c(0, 700, 1000), 3),
                    mean = 5, min = c(0, 2, 3, 0, 2, 3, 0, 2, 3), max = c(100, 8, 9, 100, 8, 9, 100, 8, 9))
  cc <- data.frame(poly = "E", metric = c("a", "b", "c"), mean = c(9, 50, 1))
  r <- nrvComparison(env, cc, c(700, 1000))
  r <- r[order(r$metric), ]
  expect_identical(r$position, c("inside", "above", "below")) ## 9 is the top of the range: inside
  expect_equal(r$lo, rep(2, 3))
  expect_equal(r$hi, rep(9, 3))
})

test_that("an all-NA NRV year does not blank the range, whichever of min, max and mean is missing", {
  env <- data.frame(poly = "E", metric = "a", time = c(700, 800, 900, 1000),
                    mean = c(5, NA, 5, 5), min = c(2, NA, NA, 3), max = c(8, NA, 9, NA))
  r <- moduleFn(".nrvRange", new.env())(env, c(700, 1000))
  expect_equal(c(r$lo, r$hi, r$mid), c(2, 9, 5))
})

test_that("stand-level labels name the leading type by sppEquiv", {
  se <- data.frame(LandR = "Pice_mar", EN_generic_full = "Black spruce")
  expect_identical(standLabel(c("lead_Pice_mar", "age_Old", "forested", "biomassTg"), se),
                   c("Black spruce leading, share of the area (%)", "Old forest, share of forested area (%)",
                     "Forested share of the area (%)", "Total tree biomass (Tg)"))
})

test_that("the stand-level figure draws reps, their median and the current condition", {
  s <- data.frame(rep = rep(1:2, each = 3), time = rep(c(0, 500, 1000), 2), class = "Old", value = 1:6)
  g <- standSeriesPanels(list(series = s, cc = data.frame(time = 0, class = "Old", value = 2),
                              period = c(700, 1000), title = "t", ylab = "y"))
  expect_s3_class(g, "ggplot")
  expect_silent(ggplot2::ggplot_build(g))
})

test_that("legend entries are the class only, on the right, in columns that keep them on the figure", {
  se <- data.frame(LandR = "Pinu_con", EN_generic_full = "Lodgepole pine")
  expect_identical(standClass(c("lead_Pinu_con", "lead_Mixed", "age_Old", "forested", "biomassTg"), se),
                   c("Lodgepole pine", "Mixed", "Old", "Forested", "Total tree biomass"))
  expect_identical(standClass(c("lead_Pinu_con", "lead_Mixed")), c("Pinu con", "Mixed")) ## no sppEquiv: one name each
  cls <- paste0("Species number ", 1:30)
  s <- expand.grid(rep = 1:2, time = c(0, 500), class = cls, stringsAsFactors = FALSE)
  s$value <- 1
  g <- standSeriesPanels(list(series = s, cc = data.frame(time = 0, class = cls, value = 1),
                              period = c(300, 500), title = "t", ylab = "y"))
  expect_identical(g$theme$legend.position, "right")
  expect_identical(g$guides$guides$colour$params$ncol, 3) ## 30 entries: 12 to a column
  expect_identical(as.character(unique(ggplot2::ggplot_build(g)$plot$data$class))[1], "Species number 1")
})

test_that("a position colours the current condition: inside green, below blue, above purple", {
  expect_identical(unname(positionColour()[c("inside", "below", "above")]), c("#1b9e77", "#2166ac", "#7b3294"))
  s <- data.frame(rep = rep(1:2, each = 3), time = rep(c(0, 500, 1000), 2), class = "Old", value = 1:6)
  g <- standSeriesPanels(list(series = s, cc = data.frame(time = 0, class = "Old", value = 2, position = "above"),
                              period = c(700, 1000), title = "t", ylab = "y"))
  b <- ggplot2::ggplot_build(g)
  pt <- b$data[[length(b$data)]]
  expect_identical(pt$fill, "#7b3294")
})

test_that("each leading type, the mixed class included, gets a figure with its position", {
  g <- expand.grid(rep = 1:5, time = c(0, 100), metric = c("lead_A", "lead_Mixed", "lead_Z", "age_Old"),
                   stringsAsFactors = FALSE)
  g$value <- ifelse(g$metric == "lead_Z", 0, 5)
  g$class <- standClass(g$metric)
  cc <- data.frame(metric = c("lead_A", "lead_Mixed", "lead_Z", "age_Old"), time = 0, value = c(1, 5, 0, 5))
  cc$class <- standClass(cc$metric)
  cmp <- data.frame(metric = cc$metric, position = c("below", "inside", "inside", "inside"))
  f <- moduleFn(".speciesFigData", new.env())(g, cc, cmp, c(50, 100), "ELF")
  expect_identical(names(f), c("lead_A", "lead_Mixed")) ## not an age class, nor a type that never leads
  expect_identical(f$lead_A$cc$position, "below")
  expect_identical(unique(f$lead_A$series$rep), 1:5)
  expect_silent(ggplot2::ggplot_build(standSeriesPanels(f$lead_A)))
})

test_that("the last stand page draws only its own panels", {
  ## 30 stand-level metrics: a page of 24 and a page of 6
  mets <- paste0("lead_S", sprintf("%02d", 1:30))
  cmp <- data.frame(poly = "SA", metric = mets, lo = 0, hi = 10, mid = 5, current = 5, position = "inside")
  pages <- moduleFn(".standOverviewPages", new.env())(cmp)
  expect_length(pages, 2L)
  expect_equal(vapply(pages, function(p) length(unique(p$panel)), integer(1)), c(24L, 6L))
  expect_false(any(pages[[2]]$panel %in% pages[[1]]$panel))
  expect_identical(sort(unlist(lapply(pages, function(p) as.character(unique(p$panel))))),
                   sort(unique(standLabel(mets))))
  b <- ggplot2::ggplot_build(moduleFn(".currentConditionOverview", new.env())(pages[[2]]))
  expect_equal(length(unique(b$layout$layout$PANEL)), 6L)
})

test_that("a comparison table is written only when asked for, beside the lm csvs", {
  withr::local_package("SpaDES.core")
  sim <- stSim()
  env <- data.frame(poly = "E", metric = "a", time = c(700, 1000), mean = 5, min = 2, max = 9)
  cc <- data.frame(poly = "E", metric = "a", mean = 9)
  write <- function(events) {
    m <- new.env()
    m$ppFiles <- character(0)
    fn <- stP(moduleFn(".writeNrvComparisonCSV", m), postprocessEvents = events)
    environment(fn)$.nrvPeriod <- function(sim) c(700, 1000)
    environment(fn)$.wantNrvComparison <- function(sim) "st" %in% events
    fn(sim, env, cc, "lm", "layerA")
    m$ppFiles
  }
  csv <- file.path(dirname(SpaDES.core::outputPath(sim)), "postprocess", "csv", "lm", "layerA", "lm_nrvComparison.csv")
  expect_length(write(c("lm", "pm", "lw", "am")), 0L) ## LandWeb's events
  expect_false(file.exists(csv))
  expect_length(write(c("lm", "st")), 1L)
  expect_true(file.exists(csv))
  expect_identical(utils::read.csv(csv)$position, "inside")
})

test_that("the comparison tables are wanted through one helper", {
  want <- function(events) stP(moduleFn(".wantNrvComparison", new.env()), postprocessEvents = events)(NULL)
  expect_false(want(c("lm", "pm", "lw", "am")))
  expect_true(want(c("lm", "st")))
  expect_true(want("ST"))
})

## a mode = "multi" sim of the module with its init scheduled but not run; the replicates' files are
## empty stand-ins for each year in `yrs`; `events` is postprocessEvents (NULL: the parameter's default)
multiSimSt <- function(events = NULL, yrs = c(700L, 1000L), interval = 300L, .env = parent.frame()) {
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
    list(mode = "multi", reps = 1L, summaryPeriod = c(700L, 1000L), summaryInterval = interval,
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

test_that("st is not run unless postprocessEvents asks", {
  withr::local_package("SpaDES.core")
  scheduled <- function(sim) {
    sim <- suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE))
    SpaDES.core::events(sim)$eventType
  }
  expect_false("st" %in% SpaDES.core::params(multiSimSt())$NRV_summary$postprocessEvents)
  expect_false("postprocess_st" %in% scheduled(multiSimSt()))
  expect_false("postprocess_st" %in% scheduled(multiSimSt(events = c("lm", "pm", "lw", "am"))))
  expect_true("postprocess_st" %in% scheduled(multiSimSt(events = c("lm", "st"))))
  expect_true("postprocess_st" %in% scheduled(multiSimSt(events = "ST")))
})

test_that("'st' is an accepted postprocess event, in any case", {
  expect_no_error(moduleFn(".checkPostprocessEvents", new.env())(c("lm", "St")))
})

test_that("with LandWeb's parameters, init schedules the events it did before", {
  withr::local_package("SpaDES.core")
  yrs <- seq(700L, 1000L, 50L)
  sim <- suppressWarnings(SpaDES.core::spades(
    multiSimSt(events = c("lm", "pm", "lw", "am"), yrs = yrs, interval = 50L),
    events = list(NRV_summary = "init"), debug = FALSE
  ))
  e <- as.data.frame(SpaDES.core::events(sim))[, c("eventTime", "eventType", "eventPriority")]
  e$eventTime <- as.numeric(e$eventTime)
  e$eventPriority <- as.numeric(e$eventPriority)
  rownames(e) <- NULL
  e <- e[order(e$eventTime, e$eventPriority), ]
  expect_equal(e, data.frame(
    eventTime = 1000,
    eventType = c("animation", "postprocess_lm", "postprocess_pm", "postprocess_lw", "postprocess", "plot"),
    eventPriority = as.numeric(SpaDES.core::.last())
  ))
})

test_that("the st figures folder is cleared only when the st event runs", {
  withr::local_package("SpaDES.core")
  sim <- stSim()
  clear <- moduleFn(".clearFigDirs", new.env())
  f <- file.path(moduleFn(".ppFigDir", new.env())(sim, "st", "SA"), "stand_age.png")
  file.create(f)
  clear(sim, "lm")
  expect_true(file.exists(f))
  clear(sim, c("lm", "st"))
  expect_false(file.exists(f))
})

test_that("the st event writes its tables and figures from snapshots on disk, and its figure folder is its own", {
  withr::local_package("SpaDES.core")
  root <- withr::local_tempdir()
  crs <- "EPSG:3400"
  rtm <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, crs = crs)
  yrs <- c(0L, 1000L)
  cd <- data.table::data.table(pixelGroup = c(1L, 2L), speciesCode = factor(c("Pice_gla", "Popu_tre")),
                               age = c(50L, 150L), B = c(2000L, 4000L))
  for (r in 1:2) {
    d <- file.path(root, "mainSim", sprintf("rep%02d", r))
    dir.create(d, recursive = TRUE)
    for (y in yrs) {
      sam <- terra::rast(rtm, vals = rep(c(50, 150), 8))
      vtm <- terra::rast(rtm, vals = rep(1:2, 8))
      levels(vtm) <- data.frame(ID = 1:2, code = c("Pice_gla", "Popu_tre"))
      terra::writeRaster(sam, file.path(d, sprintf("standAgeMap_year%04d.tif", y)))
      terra::writeRaster(vtm, file.path(d, sprintf("vegTypeMap_year%04d.tif", y)), datatype = "INT1U")
      terra::writeRaster(terra::rast(rtm, vals = rep(1:2, 8)), file.path(d, sprintf("pixelGroupMap_year%04d.tif", y)))
      qs2::qs_save(data.table::copy(cd), file.path(d, sprintf("cohortData_year%04d.qs2", y)))
    }
    file.create(file.path(d, "flammableMap_year1000.tif"))
  }
  se <- LandR::sppEquivalencies_CA
  se <- se[se[["LandR"]] %in% c("Pice_gla", "Popu_tre"), ]
  sim <- SpaDES.core::simInit(
    times = list(start = 0, end = 1000),
    params = list(NRV_summary = list(
      mode = "multi", reps = 1:2, summaryPeriod = c(0L, 1000L), summaryInterval = 1000L, timeSeriesTimes = 1:2,
      simTimes = c(0, 1000), postprocessEvents = "st", .studyAreaName = "SA", ageClassCutOffs = c(0L, 100L),
      ageClasses = c("Young", "Old")
    )),
    modules = moduleName,
    objects = list(reportingPolygons = list(A = "placeholder"), sppEquiv = se,
                   studyAreaReporting = terra::as.polygons(terra::ext(rtm), crs = crs)),
    paths = list(modulePath = modulePath, outputPath = file.path(root, "mainSim"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  )
  figDir <- file.path(root, "postprocess", "figures", "st", "SA")
  dir.create(figDir, recursive = TRUE)
  file.create(file.path(figDir, "stale.png"))
  done <- suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = c("init", "postprocess_st")),
                                              debug = FALSE))
  csv <- file.path(root, "postprocess", "csv", "st", "SA")
  st <- utils::read.csv(file.path(csv, "st.csv"))
  expect_setequal(st$rep, 1:2)
  expect_setequal(st$time, c(0, 1000))
  v <- function(m) unique(st$value[st$metric == m])
  expect_equal(v("age_Young"), 50)
  expect_equal(v("forested"), 100)
  expect_equal(v("lead_Pice_gla"), 50)
  ## 8 pixels of 2000 g/m2 and 8 of 4000, each 250 m x 250 m = 6.25 ha
  expect_equal(v("biomassTg"), (8 * 2000 + 8 * 4000) * 0.01 * 6.25 / 1e6)
  expect_true(all(file.exists(file.path(csv, c("st_CC.csv", "st_nrvComparison.csv")))))
  expect_true(all(utils::read.csv(file.path(csv, "st_nrvComparison.csv"))$position == "inside"))
  expect_false(file.exists(file.path(figDir, "stale.png")))
  expect_true(all(file.exists(file.path(figDir, c("stand_age.png", "stand_leading.png", "stand_biomass.png",
                                                  "current_condition_overview_p1.png",
                                                  file.path("species", "lead_Pice_gla.png"))))))
  expect_true(any(basename(SpaDES.core::outputs(done)$file) == "stand_biomass.png"))
})
