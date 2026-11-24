## The NRV years (`nrvWindow`): the NRV is the last fraction of `summaryPeriod`; 1 (the default) is
## all of it. The pooled lw distributions and the per-class pm pages use only those years, and the
## envelope figures shade them.

test_that(".nrvYears are the summary times in the last nrvWindow of the summary period", {
  mod <- new.env()
  nrvYears <- moduleFn(".nrvYears", mod)
  expect_identical(nrvYears(c(700, 1000), 50, 1), seq(700L, 1000L, 50L))
  expect_identical(nrvYears(c(700, 1000), 50, 0.3), c(950L, 1000L))
  expect_identical(nrvYears(c(700L, 1000L), 100L, 0.5), c(900L, 1000L))
})

test_that("reuse recomputes when the NRV years change", {
  mod <- new.env()
  aggComplete <- moduleFn(".aggComplete", mod)
  writeRepDataset <- moduleFn(".writeRepDataset", mod)
  ## the producers' output is already tidy here
  environment(writeRepDataset)$tidy_nrv_metrics <- function(x, ...) x
  environment(writeRepDataset)$write_nrv_parquet <- nrvtools::write_nrv_parquet
  root <- file.path(withr::local_tempdir(), "lw_x")
  calls <- 0L
  compute <- function(repID) {
    calls <<- calls + 1L
    data.frame(level = "landscape", class = NA_character_, id = NA_integer_, metric = "leadingProp",
               value = 1, rep = repID, time = 700L, poly = "A", metric.1 = "Pine")
  }
  build <- function(years) {
    writeRepDataset(root, "rep01", compute, reuse = TRUE, nrvYears = years)
  }

  build(c(700L, 800L))
  expect_equal(calls, 1L)
  expect_equal(readLines(file.path(root, "nrvYears.txt")), "700,800")
  expect_true(aggComplete(root, "rep01", c(700L, 800L)))
  build(c(700L, 800L)) ## same years: reused
  expect_equal(calls, 1L)
  expect_false(aggComplete(root, "rep01", 800L))
  build(800L) ## the NRV years changed: recomputed, and restamped
  expect_equal(calls, 2L)
  expect_equal(readLines(file.path(root, "nrvYears.txt")), "800")
  ## a dataset without a stamp (an earlier run's) is not reused when the NRV years are given
  unlink(file.path(root, "nrvYears.txt"))
  expect_false(aggComplete(root, "rep01", 800L))
  ## no NRV years given: the stamp is not looked at
  expect_true(aggComplete(root, "rep01"))
})

test_that("pm by-class figures pool the NRV years", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  aggRoot <- moduleFn(".nrvAggRoot", mod)
  summarizeRepDataset <- moduleFn(".summarizeRepDataset", mod)
  refCode <- LandWebUtils::refCodeFor("pm", "FMA")
  unlink(c(aggRoot(sim, refCode), aggRoot(sim, paste0(refCode, "_CC"))), recursive = TRUE)
  one <- function(r, time) {
    data.frame(level = "class", class = c("Decid", "Pine"), id = NA_integer_, metric = "ca",
               value = c(10, 20) * r, rep = r, time = time, poly = "A")
  }
  for (r in 1:2) {
    nrvtools::write_nrv_parquet(do.call(rbind, lapply(c(700L, 800L, 900L, 1000L), function(t) one(r, t))),
                                aggRoot(sim, refCode), sprintf("rep%02d", r))
  }
  nrvtools::write_nrv_parquet(one(1, 0L), aggRoot(sim, paste0(refCode, "_CC")), "CC")
  mod[[refCode]] <- summarizeRepDataset(aggRoot(sim, refCode), perReplicate = TRUE)
  mod[[paste0(refCode, "_CC")]] <- summarizeRepDataset(aggRoot(sim, paste0(refCode, "_CC")), perReplicate = TRUE)
  byClassFigTasks <- moduleFn(".byClassFigTasks", mod)

  all <- byClassFigTasks(sim, "pm", "FMA")[[1L]]
  nrv <- byClassFigTasks(sim, "pm", "FMA", nrvPeriod = c(900, 1000))[[1L]]

  expect_setequal(all$df$time, c(700, 800, 900, 1000))
  expect_setequal(nrv$df$time, c(900, 1000))
  expect_equal(nrow(nrv$df), 2L * 2L * 2L) ## replicates x classes x NRV years
  expect_equal(nrow(nrv$cc), 2L) ## the current condition is one snapshot, whatever the NRV years
})

test_that("shading reaches the render workers", {
  ## render_one runs detached from the module namespace, so every module helper it calls must be
  ## re-assigned into its environment by .renderTasks(). A package attached on the search path
  ## would hide a missing one (the detached environment sits under globalenv()), so the test runs
  ## with it detached, as in a simulation.
  mod <- new.env()
  renderTasks <- moduleFn(".renderTasks", mod)
  environment(renderTasks)$.plotWorkers <- function(sim, nTasks) 1L
  g <- expand.grid(time = c(100, 200, 300), class = "Pine", stringsAsFactors = FALSE)
  df <- data.frame(time = g$time, poly = "ELF", class = g$class, metric = "area_mn", metric.1 = NA_character_,
                   level = "class", n_reps = 3L, mean = 10, min = 8, max = 12, q25 = 9, median = 10, q75 = 11)
  cc <- data.frame(poly = "ELF", class = "Pine", metric = "area_mn", metric.1 = NA_character_, mean = 11)
  file <- file.path(withr::local_tempdir(), "env.png")
  task <- list(plotter = "envelope", df = df, type = "ribbon", facet = c("class", "metric.1"),
               ylab = "y", title = "t", page = NULL, nrvPeriod = c(200, 300), cc = cc,
               file = file, width = 5, height = 4)
  fn <- get0(".renderTasks", mode = "function")
  if (!is.null(fn) && isNamespace(topenv(environment(fn)))) {
    pkgEnv <- paste0("package:", getNamespaceName(topenv(environment(fn))))
    if (pkgEnv %in% search()) {
      attached <- as.environment(pkgEnv)
      pos <- match(pkgEnv, search())
      detach(pkgEnv, character.only = TRUE)
      withr::defer(attach(attached, pos = pos, name = pkgEnv, warn.conflicts = FALSE))
    }
  }

  expect_no_warning(out <- renderTasks(NULL, list(task)))
  expect_equal(out, file)
  expect_true(file.exists(file))
  ## the shading is drawn: the same figure without NRV years differs
  plain <- task
  plain$nrvPeriod <- NULL
  plain$file <- file.path(dirname(file), "plain.png")
  renderTasks(NULL, list(plain))
  expect_false(unname(tools::md5sum(file) == tools::md5sum(plain$file)))
})

test_that(".addNrvShading shades the NRV years on both time axes and keeps the caption", {
  mod <- new.env()
  addNrvShading <- moduleFn(".addNrvShading", mod)
  g <- expand.grid(time = c(100, 200, 300), class = "Pine", stringsAsFactors = FALSE)
  df <- data.frame(time = g$time, poly = "ELF", class = g$class, metric = "area_mn", metric.1 = NA_character_,
                   level = "class", n_reps = 3L, mean = 10, min = 8, max = 12, q25 = 9, median = 10, q75 = 11)
  for (type in c("ribbon", "boxplot")) {
    gg <- addNrvShading(nrvtools::plot_nrv_envelope(df, type = type, facet = "class"), c(200, 300))
    rect <- ggplot2::ggplot_build(gg)$data[[1]] ## drawn first, behind the data
    expect_match(gg$labels$caption, "shaded: NRV years")
    expect_equal(unique(as.numeric(c(rect$xmin, rect$xmax))), if (type == "boxplot") c(1.5, 3.5) else c(200, 300))
  }
  ## no NRV years given: the figure is returned as it is
  gg <- nrvtools::plot_nrv_envelope(df, type = "ribbon", facet = "class")
  expect_identical(addNrvShading(gg, NULL), gg)
})

test_that("init stops when no summary time falls in the NRV years", {
  pgm <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, vals = 1)
  names(pgm) <- "pixelGroup"
  sa <- terra::as.polygons(terra::ext(pgm))
  terra::crs(sa) <- terra::crs(pgm)
  objs <- list(
    cohortData = data.table::data.table(pixelGroup = 1L, speciesCode = "Pice_mar", age = 30L, B = 1000L),
    pixelGroupMap = pgm, studyAreaReporting = sa, speciesLayers = terra::rast(pgm, vals = 1),
    sppColorVect = c(Pice_mar = "black", Mixed = "grey"),
    sppEquiv = data.table::data.table(LandR = "Pice_mar", Type = "Conifer", EN_generic_full = "Black spruce")
  )
  root <- withr::local_tempdir()
  sim <- SpaDES.core::simInit(
    times = list(start = 0, end = 1000),
    params = list(NRV_summary = list(mode = "single", summaryPeriod = c(700, 1000), summaryInterval = 200L,
                                     nrvWindow = 0.1, timeSeriesTimes = c(0, 1000), .plots = NA)),
    modules = moduleName, objects = objs,
    paths = list(modulePath = modulePath, outputPath = file.path(root, "out"),
                 cachePath = file.path(root, "cache"), inputPath = file.path(root, "in"))
  )
  expect_error(
    suppressWarnings(SpaDES.core::spades(sim, events = list(NRV_summary = "init"), debug = FALSE)),
    "no summary time falls in the NRV years.*nrvWindow 0.1"
  )
})
