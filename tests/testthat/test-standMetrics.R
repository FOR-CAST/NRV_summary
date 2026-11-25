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

test_that("snapshot files pair by rep folder and year, the parquet cohorts too", {
  f <- function(kind) sprintf("o/rep%d/%s_year%04d.tif", c(2, 1, 1), kind, c(2120, 2020, 2120))
  cd <- sprintf("o/rep%d/cohortData/year=%04d/part-0.parquet", c(1, 2, 1), c(2120, 2120, 2020))
  s <- moduleFn(".snapshotTable", new.env())(f("standAgeMap"), f("vegTypeMap"), f("pixelGroupMap"), cd)
  expect_identical(nrow(s), 3L)
  i <- s$rep == "rep2"
  expect_identical(s$time[i], 2120)
  expect_identical(s$cd[i], "o/rep2/cohortData/year=2120/part-0.parquet")
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

test_that("a comparison table is written only when asked for, beside the lm csvs, landscape level only", {
  withr::local_package("SpaDES.core")
  sim <- stSim()
  ## a class-level row of the same metric, with another range, must not enter the table
  env <- data.frame(poly = "E", level = rep(c("landscape", "class"), each = 2), class = rep(c(NA, "X"), each = 2),
                    metric = "a", time = c(700, 1000), mean = 5, min = rep(c(2, 0), each = 2),
                    max = rep(c(9, 100), each = 2))
  cc <- data.frame(poly = "E", level = c("landscape", "class"), class = c(NA, "X"), metric = "a", mean = c(9, 50))
  recorded <- character(0)
  write <- function(events, layer = "layerA") {
    m <- new.env()
    m$ppFiles <- character(0)
    fn <- stP(moduleFn(".writeNrvComparisonCSV", m), postprocessEvents = events)
    environment(fn)$.nrvPeriod <- function(sim) c(700, 1000)
    environment(fn)$.wantNrvComparison <- function(sim) any(c("st", "cc") %in% events)
    environment(fn)$.recordComparison <- function(kind, cmp, layer) recorded <<- c(recorded, paste(kind, layer))
    fn(sim, env, cc, "lm", layer)
    m$ppFiles
  }
  csvDir <- function(layer) file.path(dirname(SpaDES.core::outputPath(sim)), "postprocess", "csv", "lm", layer)
  csv <- file.path(csvDir("layerA"), "lm_nrvComparison.csv")
  expect_length(write(c("lm", "pm", "lw", "am")), 0L) ## LandWeb's events
  expect_false(file.exists(csv))
  expect_false(dir.exists(csvDir("layerA"))) ## not even the folder
  expect_length(recorded, 0L) ## nothing is kept for the combined figure either
  expect_length(write(c("lm", "st")), 1L)
  expect_identical(recorded, "lm layerA") ## the figure of every group gets it
  expect_length(write(c("lm", "cc")), 1L) ## cc asks for the table too
  expect_true(file.exists(csv))
  r <- utils::read.csv(csv)
  expect_identical(r$position, "inside")
  expect_equal(c(nrow(r), r$lo, r$hi), c(1, 2, 9)) ## the landscape row's range, not the class row's
  ## a table of an earlier run with st does not stay beside the csvs of one without it
  expect_length(write(c("lm", "pm")), 0L)
  expect_false(file.exists(csv))
})

test_that("the comparison tables are wanted through one helper", {
  want <- function(events) stP(moduleFn(".wantNrvComparison", new.env()), postprocessEvents = events)(NULL)
  expect_false(want(c("lm", "pm", "lw", "am")))
  expect_true(want(c("lm", "st")))
  expect_true(want("ST"))
})

## a mode = "multi" sim of the module with its init scheduled but not run; the replicates' files are
## empty stand-ins for each year in `yrs`; `events` is postprocessEvents (NULL: the parameter's default)
multiSimSt <- function(events = NULL, yrs = c(700L, 1000L), interval = 300L, params = list(),
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
    list(mode = "multi", reps = 1L, summaryPeriod = c(700L, 1000L), summaryInterval = interval,
         timeSeriesTimes = 701:702, simTimes = c(0, 1000)),
    if (!is.null(events)) list(postprocessEvents = events),
    params
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

## two reps of a 4 x 4 landscape (two species, two stand ages) at years 0 and 1000 on disk, and the sim of
## the module over them with `events` as postprocessEvents; the reporting polygon is the whole landscape.
## `cohorts`: the format of the saved cohortData, "qs2", "parquet", "both" (the .qs2 files then hold other cohorts,
## with ten times the biomass and the species swapped, so the results show which were read) or "stale" (the .qs2
## files hold the cohorts and the parquet ones the others). `outputsDF`: the files are given as an outputsDF, which
## with "stale" does not list the parquet ones. `reporting`: the species are reported by group (Spruce, Aspen).
## `summaryPeriod`: the module's parameter (the default makes year 0 a summary year, which InitMulti() then checks).
stDiskSim <- function(events, subregion = "all", studyArea = "SA", cohorts = "qs2", outputsDF = FALSE,
                      reporting = FALSE, summaryPeriod = c(0L, 1000L), .env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .env)
  crs <- "EPSG:3400"
  rtm <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, crs = crs)
  cd <- data.table::data.table(pixelGroup = c(1L, 2L), speciesCode = factor(c("Pice_gla", "Popu_tre")),
                               age = c(50L, 150L), B = c(2000L, 4000L))
  other <- data.table::copy(cd)
  data.table::set(other, j = "B", value = other$B * 10L)
  data.table::set(other, j = "speciesCode", value = rev(other$speciesCode))
  qsCd <- if (cohorts == "both") other else cd
  pqCd <- if (cohorts == "stale") other else cd
  for (r in 1:2) {
    d <- file.path(root, "mainSim", sprintf("rep%02d", r))
    dir.create(d, recursive = TRUE)
    for (y in c(0L, 1000L)) {
      sam <- terra::rast(rtm, vals = rep(c(50, 150), 8))
      vtm <- terra::rast(rtm, vals = rep(1:2, 8))
      levels(vtm) <- data.frame(ID = 1:2, code = c("Pice_gla", "Popu_tre"))
      terra::writeRaster(sam, file.path(d, sprintf("standAgeMap_year%04d.tif", y)))
      terra::writeRaster(vtm, file.path(d, sprintf("vegTypeMap_year%04d.tif", y)), datatype = "INT1U")
      terra::writeRaster(terra::rast(rtm, vals = rep(1:2, 8)), file.path(d, sprintf("pixelGroupMap_year%04d.tif", y)))
      if (cohorts %in% c("qs2", "both", "stale")) {
        qs2::qs_save(data.table::copy(qsCd), file.path(d, sprintf("cohortData_year%04d.qs2", y)))
      }
      if (cohorts %in% c("parquet", "both", "stale")) {
        nrvtools::write_cohort_parquet(pqCd, file.path(d, "cohortData"), y)
      }
    }
    file.create(file.path(d, "flammableMap_year1000.tif"))
  }
  se <- LandR::sppEquivalencies_CA
  se <- data.table::copy(se[se[["LandR"]] %in% c("Pice_gla", "Popu_tre"), ])
  data.table::set(se, j = "Group", value = ifelse(se[["LandR"]] == "Pice_gla", "Spruce", "Aspen"))
  listed <- list.files(file.path(root, "mainSim"), recursive = TRUE, full.names = TRUE)
  if (cohorts == "stale") listed <- listed[basename(listed) != "part-0.parquet"]
  poly <- terra::as.polygons(terra::ext(rtm), crs = crs)
  poly$Name <- subregion
  list(root = root, sim = SpaDES.core::simInit(
    times = list(start = 0, end = 1000),
    params = list(NRV_summary = list(
      mode = "multi", reps = 1:2, summaryPeriod = summaryPeriod, summaryInterval = 1000L, timeSeriesTimes = 1:2,
      simTimes = c(0, 1000), postprocessEvents = events, .studyAreaName = studyArea, ageClassCutOffs = c(0L, 100L),
      ageClasses = c("Young", "Old"), sppEquivColReporting = if (reporting) "Group" else NA_character_
    )),
    modules = moduleName,
    objects = c(list(reportingPolygons = list(A = sf::st_as_sf(poly)), sppEquiv = se, studyAreaReporting = poly),
                if (outputsDF) list(outputsDF = data.frame(file = listed))),
    paths = list(modulePath = modulePath, outputPath = file.path(root, "mainSim"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  ))
}

test_that("the st event writes its tables and figures from snapshots on disk, and its figure folder is its own", {
  withr::local_package("SpaDES.core")
  fx <- stDiskSim(c("st", "stability"))
  root <- fx$root
  figDir <- file.path(root, "postprocess", "figures", "st", "SA")
  dir.create(figDir, recursive = TRUE)
  file.create(file.path(figDir, "stale.png"))
  done <- suppressWarnings(SpaDES.core::spades(
    fx$sim, debug = FALSE,
    events = list(NRV_summary = c("init", "postprocess_st", "postprocess", "plot", "postprocess_stability"))
  ))
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
  pngs <- file.path(figDir, c("stand_age.png", "stand_leading.png", "stand_biomass.png",
                              "current_condition_overview_p1.png", file.path("species", "lead_Pice_gla.png")))
  expect_true(all(file.exists(pngs))) ## still there after the plot event
  out <- SpaDES.core::outputs(done)$file
  expect_true(all(c("stand_biomass.png", "st.csv", "st_CC.csv", "st_nrvComparison.csv") %in% basename(out)))
  ## with stability in the same run, the stand series are tested too
  stab <- utils::read.csv(file.path(root, "postprocess", "csv", "stability", "all", "stability.csv"))
  expect_true("st" %in% stab$kind)
  expect_true(all(c("age_Young", "biomassTg") %in% stab$metric[stab$kind == "st"]))
})

test_that("with LandWeb's events the lm event writes no comparison table, with st it does", {
  withr::local_package("SpaDES.core")
  run <- function(events) {
    fx <- stDiskSim(events)
    suppressWarnings(suppressMessages(SpaDES.core::spades(
      fx$sim, debug = FALSE, events = list(NRV_summary = c("init", "postprocess_lm"))
    )))
    list.files(file.path(fx$root, "postprocess", "csv", "lm"), recursive = TRUE)
  }
  without <- run(c("lm", "pm", "lw", "am"))
  expect_true(length(without) > 0L) ## the lm event ran
  expect_false(any(grepl("nrvComparison", without)))
  with <- run(c("lm", "st"))
  expect_identical(setdiff(with, without), file.path("A", "lm_nrvComparison.csv"))
  withCc <- run(c("lm", "cc")) ## cc asks for the tables too
  expect_identical(setdiff(withCc, without), file.path("A", "lm_nrvComparison.csv"))
})

test_that("a code shared by several sppEquiv rows is named by all of them", {
  se <- data.frame(LandR = c("Pice_gla", "Pice_gla", "Betu_pap"),
                   EN_generic_full = c("White spruce", "Engelmann spruce", "Paper birch"))
  expect_identical(standClass(c("lead_Pice_gla", "lead_Betu_pap", "lead_Mixed"), se),
                   c("White spruce / Engelmann spruce", "Paper birch", "Mixed"))
  expect_identical(moduleFn(".annualSeriesLabel", new.env())("lead_Pice_gla", se),
                   "Leading: White spruce / Engelmann spruce")
})

test_that("an empty figure group is skipped, not an error", {
  expect_identical(moduleFn(".standGroups", new.env())(c(lead_A = 0.2, lead_B = 0.1), c("Young", "Old"))$leading$m,
                   character(0))
  expect_identical(names(Filter(function(g) length(g$m) > 0L,
                                moduleFn(".standGroups", new.env())(c(lead_A = 0.2), c("Young")))),
                   c("age", "biomass"))
})

test_that("the stand series join the stability check only when both st and stability are requested", {
  m <- new.env()
  m$rptPolyNames <- character(0)
  m$standSeries <- data.frame(metric = rep(c("age_Old", "lead_A"), each = 6), rep = rep(1:2, each = 3),
                              time = c(0, 500, 1000), value = 1:12)
  series <- function(events) {
    fn <- stP(moduleFn(".stabilitySeries", m), postprocessEvents = events)
    environment(fn)$.stabilityBurn <- function(sim) NULL
    environment(fn)$.studyAreaLabel <- function(sim) "SA"
    suppressMessages(fn(NULL))
  }
  r <- series(c("st", "stability"))
  expect_setequal(r$kind, "st")
  expect_setequal(r$metric, c("age_Old", "lead_A"))
  expect_setequal(c(r$layer, r$poly), "SA")
  expect_setequal(r$level, "landscape")
  expect_false("st" %in% series("stability")$kind) ## a left-over series of an earlier st event is not tested
  expect_identical(moduleFn(".stabilityTitle", new.env())(c("st", "st"), c("age_Old", "forested"), "landscape"),
                   standLabel(c("age_Old", "forested")))
})

## ---- the opt-in `cc` event: scheduling, the comparison tables and a run (the rest is in test-current-condition-all.R)

test_that("cc is not drawn unless postprocessEvents asks", {
  withr::local_package("SpaDES.core")
  scheduled <- function(sim) {
    sim <- suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE))
    SpaDES.core::events(sim)$eventType
  }
  expect_false("postprocess_cc" %in% scheduled(multiSimSt()))
  expect_false("postprocess_cc" %in% scheduled(multiSimSt(events = c("lm", "pm", "lw", "am", "st"))))
  expect_true("postprocess_cc" %in% scheduled(multiSimSt(events = c("lm", "cc"))))
  expect_true("postprocess_cc" %in% scheduled(multiSimSt(events = c("lm", "CC"))))
  expect_no_error(moduleFn(".checkPostprocessEvents", new.env())(c("lm", "cc")))
})

test_that("cc runs after the events that record the comparisons and before the plot event", {
  withr::local_package("SpaDES.core")
  sim <- suppressWarnings(SpaDES.core::spades(
    multiSimSt(events = c("lm", "st", "autocorr", "cc", "stability")), events = list(NRV_summary = "init"), debug = FALSE
  ))
  e <- as.data.frame(SpaDES.core::events(sim))
  e <- e[order(as.numeric(e$eventTime), as.numeric(e$eventPriority)), ]
  ev <- e$eventType[e$eventType %in% c("postprocess_lm", "postprocess_st", "postprocess_autocorr", "postprocess_cc",
                                       "plot", "postprocess_stability")]
  expect_true(all(c("postprocess_autocorr", "postprocess_stability") %in% ev))
  expect_identical(ev[1:2], c("postprocess_lm", "postprocess_st"))
  expect_lt(match("postprocess_st", ev), match("postprocess_cc", ev))
  expect_lt(match("postprocess_cc", ev), match("plot", ev))
  expect_lt(match("plot", ev), match("postprocess_stability", ev))
})

test_that("cc writes the comparison tables", {
  want <- function(events) stP(moduleFn(".wantNrvComparison", new.env()), postprocessEvents = events)(NULL)
  expect_false(want(c("lm", "pm", "lw", "am")))
  expect_true(want(c("lm", "cc")))
  expect_true(want("CC"))
})

test_that("a run of st and cc draws the figure in the study area's layer, and the plot event keeps it", {
  withr::local_package("SpaDES.core")
  fx <- stDiskSim(c("st", "cc"))
  all <- file.path(fx$root, "postprocess", "figures", "all")
  dir.create(file.path(all, "SA"), recursive = TRUE)
  stale <- file.path(all, "SA", "stale.png")
  file.create(stale)
  done <- suppressWarnings(suppressMessages(SpaDES.core::spades(
    fx$sim, debug = FALSE,
    events = list(NRV_summary = c("init", "postprocess_st", "postprocess_cc", "postprocess", "plot"))
  )))
  expect_false(file.exists(stale))
  f <- file.path(all, "SA", "SA.png") ## the unit is the study area, in its own layer
  expect_true(file.exists(f))
  expect_true(normalizePath(f) %in% normalizePath(SpaDES.core::outputs(done)$file))
})

test_that("a run of lm, st and cc puts the stand rows on the figure of the unit named like the study area", {
  withr::local_package("SpaDES.core")
  run <- function(subregion) {
    fx <- stDiskSim(c("lm", "st", "cc"), subregion = subregion)
    suppressWarnings(suppressMessages(SpaDES.core::spades(
      fx$sim, debug = FALSE, events = list(NRV_summary = c("init", "postprocess_lm", "postprocess_st", "postprocess_cc"))
    )))
    list.files(file.path(fx$root, "postprocess", "figures", "all"), recursive = TRUE)
  }
  expect_identical(run("SA"), file.path("A", "SA.png")) ## one figure, the layer's unit with the stand rows
  expect_identical(sort(run("all")), sort(c(file.path("A", "all.png"), file.path("SA", "SA.png")))) ## LandWeb's case
})

## ---- the P5 final-review fix wave ----

test_that("a leading type found only in the current condition, or above the threshold only there, is drawn", {
  rows <- function(metric, rep, value) data.frame(metric = metric, rep = rep, time = c(0, 100), value = value)
  series <- rbind(rows("lead_A", 1, c(60, 60)), rows("lead_B", 1, c(0.3, 0.2)), rows("lead_C", 1, c(0, 0)),
                  rows("lead_D", 1, c(0, 0)))
  cc <- data.frame(metric = c("lead_A", "lead_B", "lead_C", "lead_D"), time = 0, value = c(50, 5, 4, 0.1))
  lm <- moduleFn(".leadMax", new.env())(series, cc)
  expect_equal(lm[["lead_B"]], 5)
  expect_equal(lm[["lead_C"]], 4) ## present only in the current condition
  m <- moduleFn(".standGroups", new.env())(lm, c("Young", "Old"))$leading$m
  expect_setequal(m, c("lead_A", "lead_B", "lead_C")) ## lead_D leads 0.1% at most
  ## the species figures too: a type that leads nowhere in the replicates but leads today
  series$class <- standClass(series$metric)
  cc$class <- standClass(cc$metric)
  cmp <- data.frame(metric = cc$metric, position = "above")
  f <- moduleFn(".speciesFigData", new.env())(series, cc, cmp, c(0, 100), "SA")
  expect_setequal(names(f), c("lead_A", "lead_B", "lead_C", "lead_D"))
})

test_that("the age classes are ordered young to old in the stand figure legend", {
  metric <- c("age_Old", "age_Young", "age_Mature", "lead_Pice_gla", "age_Immature", "forested", "biomassTg")
  cls <- standClass(metric)
  f <- moduleFn(".standClassFactor", new.env())(metric, cls, c("Young", "Immature", "Mature", "Old"))
  expect_s3_class(f, "factor")
  expect_identical(as.character(f), cls)
  expect_identical(levels(f)[1:6], c("Total tree biomass", "Young", "Immature", "Mature", "Old", "Forested"))
  s <- data.frame(rep = 1L, time = rep(c(0, 100), 4), class = f[rep(c(2, 5, 3, 1), each = 2)], value = 1)
  g <- standSeriesPanels(list(series = s, cc = data.frame(time = 0, class = f[2], value = 1),
                              period = c(0, 100), title = "t", ylab = "y"))
  expect_identical(ggplot2::ggplot_build(g)$plot$scales$get_scales("colour")$get_limits(),
                   c("Young", "Immature", "Mature", "Old"))
})

test_that("with ten or more age classes the forested share still follows them", {
  ac <- paste0("C", sprintf("%02d", 1:12))
  rank <- function(m) moduleFn(".standRank", new.env())(m, rep("stand", length(m)), ac)
  r <- rank(c("biomassTg", paste0("age_", ac), "forested", "lead_X"))
  expect_true(all(diff(r[1:13]) > 0)) ## biomass, then the 12 age classes in order
  expect_true(rank("forested") > rank("age_C12"))
  expect_true(rank("lead_X") > rank("forested"))
})

test_that("a missing study-area name is 'unit' everywhere the st event uses it", {
  label <- function(name) stP(moduleFn(".studyAreaLabel", new.env()), .studyAreaName = name)(NULL)
  expect_identical(label("SA"), "SA")
  for (bad in list(NA_character_, "", NULL)) expect_identical(label(bad), "unit")
})

test_that("the st event with no study-area name writes under 'unit'", {
  withr::local_package("SpaDES.core")
  fx <- stDiskSim(c("st", "stability"), studyArea = NA_character_)
  suppressWarnings(SpaDES.core::spades(
    fx$sim, debug = FALSE,
    events = list(NRV_summary = c("init", "postprocess_st", "postprocess_stability"))
  ))
  pp <- file.path(fx$root, "postprocess")
  expect_true(file.exists(file.path(pp, "csv", "st", "unit", "st.csv")))
  expect_true(file.exists(file.path(pp, "figures", "st", "unit", "stand_age.png")))
  expect_false(dir.exists(file.path(pp, "csv", "st", "NA")))
  stab <- utils::read.csv(file.path(pp, "csv", "stability", "all", "stability.csv"))
  expect_setequal(stab$layer[stab$kind == "st"], "unit")
})

test_that("a stale comparison table goes when the table is wanted but there is nothing to write", {
  withr::local_package("SpaDES.core")
  sim <- stSim()
  f <- file.path(dirname(SpaDES.core::outputPath(sim)), "postprocess", "csv", "lm", "L", "lm_nrvComparison.csv")
  dir.create(dirname(f), recursive = TRUE)
  fn <- stP(moduleFn(".writeNrvComparisonCSV", new.env()), postprocessEvents = c("lm", "st"))
  environment(fn)$.wantNrvComparison <- function(sim) TRUE
  for (empty in list(NULL, data.frame())) {
    file.create(f)
    fn(sim, empty, data.frame(poly = "E", metric = "a", mean = 1), "lm", "L")
    expect_false(file.exists(f))
  }
  env <- data.frame(poly = "E", level = "class", metric = "a", time = 700, mean = 5, min = 2, max = 9)
  file.create(f)
  fn(sim, env, data.frame(poly = "E", level = "class", metric = "a", mean = 1), "lm", "L") ## no landscape rows
  expect_false(file.exists(f))
})

test_that("recording a comparison returns nothing, for no rows as for some", {
  m <- new.env()
  m$nrvComparisons <- list()
  rec <- moduleFn(".recordComparison", m)
  cmp <- data.frame(poly = "E", metric = "ai", lo = 0, hi = 1, mid = 0.5, current = 0.5, position = "inside")
  expect_null(rec("lm", cmp, "L"))
  expect_null(rec("lm", cmp[0, ], "L"))
  expect_named(m$nrvComparisons, "lm L")
})

test_that("init stops on an envelopeStyle that is neither ribbon nor lines", {
  withr::local_package("SpaDES.core")
  sim <- multiSimSt(params = list(envelopeStyle = "dots"))
  expect_error(suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)),
               "envelopeStyle must be")
})

test_that("st gives the same tables from .qs2 and parquet cohorts", {
  withr::local_package("SpaDES.core")
  tables <- function(cohorts) {
    fx <- stDiskSim("st", cohorts = cohorts)
    suppressWarnings(SpaDES.core::spades(fx$sim, debug = FALSE,
                                         events = list(NRV_summary = c("init", "postprocess_st"))))
    csv <- file.path(fx$root, "postprocess", "csv", "st", "SA")
    lapply(c(st = "st.csv", cc = "st_CC.csv", cmp = "st_nrvComparison.csv"), function(f) readLines(file.path(csv, f)))
  }
  qs <- tables("qs2")
  expect_true(all(lengths(qs) > 1L))
  expect_identical(tables("parquet"), qs)
  expect_identical(tables("both"), qs) ## the parquet cohorts are read, at year 0 as at the summary times
})

test_that("given an outputsDF, every site reads the .qs2 cohorts it lists, not stale parquet ones beside them", {
  ## the reporting-group maps (year 0 and the summary years), the st snapshots and the current-condition row
  withr::local_package("SpaDES.core")
  run <- function(cohorts) {
    fx <- stDiskSim("st", cohorts = cohorts, outputsDF = TRUE, reporting = TRUE)
    suppressWarnings(SpaDES.core::spades(fx$sim, debug = FALSE,
                                         events = list(NRV_summary = c("init", "postprocess_st"))))
    csv <- file.path(fx$root, "postprocess", "csv", "st", "SA")
    mapDir <- file.path(fx$root, "postprocess", "_reportingMaps")
    maps <- sort(list.files(mapDir, "\\.tif$", recursive = TRUE))
    list(st = readLines(file.path(csv, "st.csv")), cc = readLines(file.path(csv, "st_CC.csv")), maps = maps,
         values = lapply(maps, function(f) terra::values(terra::rast(file.path(mapDir, f)), mat = FALSE)))
  }
  qs <- run("qs2")
  expect_true(all(c(file.path("rep01", "vegTypeMap_year0000.tif"), file.path("rep02", "vegTypeMap_year1000.tif"))
                  %in% qs$maps))
  expect_identical(run("stale"), qs)
})

test_that("st stops, naming both formats' paths, when the year-0 cohort file is missing", {
  ## summaryPeriod 1000 only: InitMulti() does not check year 0, so the event is the one to stop
  withr::local_package("SpaDES.core")
  fx <- stDiskSim("st", summaryPeriod = c(1000L, 1000L))
  rep1 <- file.path(fx$root, "mainSim", "rep01")
  expect_true(file.remove(file.path(rep1, "cohortData_year0000.qs2")))
  expect_error(
    suppressWarnings(SpaDES.core::spades(fx$sim, debug = FALSE,
                                         events = list(NRV_summary = c("init", "postprocess_st")))),
    paste0("cannot take the stand-level measures of the current condition; missing:\n",
           file.path(rep1, "cohortData_year0000.qs2"), " or ",
           file.path(rep1, "cohortData", "year=0000", "part-0.parquet")),
    fixed = TRUE
  )
})

test_that("the bc year-0 check stops, naming both formats' paths, when the year-0 cohort file is missing", {
  ## the bc event needs the BC seral inputs (landcover and seral tables), so the two steps it takes are run on their
  ## own: .repYear0Files() chooses the file and .checkCohortFile() stops
  year0 <- moduleFn(".repYear0Files", new.env())
  check <- moduleFn(".checkCohortFile", new.env())
  rep1 <- file.path(withr::local_tempdir(), "rep01")
  dir.create(rep1)
  f <- year0(rep1, file.path(rep1, "cohortData_year0700.qs2"))$cohortData
  expect_error(
    check(f, "cannot make the year-0 seral stage map"),
    paste0("cannot make the year-0 seral stage map; missing:\n", file.path(rep1, "cohortData_year0000.qs2"), " or ",
           file.path(rep1, "cohortData", "year=0000", "part-0.parquet")),
    fixed = TRUE
  )
  ## and it does not stop when either format is there
  file.create(file.path(rep1, "cohortData_year0000.qs2"))
  expect_identical(check(year0(rep1)$cohortData, "x"), file.path(rep1, "cohortData_year0000.qs2"))
})
