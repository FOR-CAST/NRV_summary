## One figure with every NRV metric of every group (the opt-in `cc` event), each metric's NRV range rescaled to
## 0-1. Ported from the PE fork (PE#21), with a figure per reporting layer and unit, and file names that cannot
## collide. Opt-in: with LandWeb's parameters a run schedules no `cc` event and writes no comparison table.

ccFn <- function(name, ...) moduleFn(name, new.env())(...)

cmpRows <- function(poly = "ELF", layer = "L1") {
  data.frame(
    layer = layer, poly = poly, metric = c("ai", "area_mn", "age_Old", "lead_A", "lead_Zero", "forested"),
    group = c("landscape", "landscape", "stand", "stand", "stand", "stand"),
    label = c("Aggregation index", "Mean patch area", "Old forest", "A leading", "Zero leading", "Forested"),
    lo = c(10, 100, 40, 1, 0, 50), hi = c(20, 200, 60, 3, 0, 50), mid = c(15, 150, 50, 2, 0, 50),
    current = c(15, 50, 90, 2, 0, 50), position = c("inside", "below", "above", "inside", "inside", "inside")
  )
}

test_that("a value is rescaled so the NRV range is 0 to 1", {
  expect_equal(ccFn(".nrvScale", c(10, 15, 20, 5, 30), 10, 20), c(0, 0.5, 1, -0.5, 2))
})

test_that("a range of one value gets a fixed small width around it", {
  x <- ccFn(".nrvScale", c(49.5, 50, 50.5), 50, 50)
  expect_equal(x[2], 0.5) ## the value is the middle of the bar
  expect_equal(ccFn(".nrvScale", c(0, 1), 0, 0), c(0.5, 100.5))
  expect_true(all(is.finite(x)))
})

test_that("rows of every group appear, and a leading type that is 0 everywhere is dropped", {
  d <- ccFn(".currentConditionAllData", cmpRows())
  expect_setequal(d$metric, c("ai", "area_mn", "age_Old", "lead_A", "forested"))
  expect_setequal(d$group, c("landscape", "stand"))
  expect_false("lead_Zero" %in% d$metric)
  expect_equal(d$xCurrent[d$metric == "age_Old"], 2.5) ## (90 - 40) / 20 ...
})

test_that("far-out points are clipped to the axis limit and keep their true value", {
  d <- ccFn(".currentConditionAllData", cmpRows())
  lim <- ccFn(".ccAllLimits")
  expect_true(all(d$xPlot >= lim[1] & d$xPlot <= lim[2]))
  expect_equal(d$xPlot[d$metric == "age_Old"], lim[2])
  expect_true(d$clipped[d$metric == "age_Old"])
  expect_false(d$clipped[d$metric == "ai"])
})

test_that("the figure is one panel, with the metrics as rows coloured by position", {
  gg <- ccFn(".currentConditionAll", ccFn(".currentConditionAllData", cmpRows()))
  b <- ggplot2::ggplot_build(gg)
  expect_s3_class(gg$facet, "FacetNull")
  expect_equal(length(unique(b$layout$layout$PANEL)), 1L)
  expect_identical(gg$labels$x, ccFn(".ccAllXlab"))
  expect_true(all(unname(ccFn(".positionColour")) %in% unlist(lapply(b$data, function(l) l$colour)))) ## inside, below, above
})

test_that("each group's comparison is recorded with its layer and label, for the combined figure", {
  mod <- new.env()
  record <- moduleFn(".recordComparison", mod)
  record("lm", data.frame(poly = "ELF", metric = "ai", lo = 1, hi = 2, mid = 1.5, current = 3, position = "above"), "L1")
  record("st", data.frame(poly = "ELF", metric = "age_Old", lo = 1, hi = 2, mid = 1.5, current = 1, position = "inside",
                          label = "Old forest"), "SA")
  d <- do.call(rbind, mod$nrvComparisons)
  expect_setequal(d$group, c("landscape", "stand"))
  expect_equal(d$layer[d$metric == "ai"], "L1")
  expect_equal(d$layer[d$metric == "age_Old"], "SA")
  expect_equal(d$label[d$metric == "age_Old"], "Old forest")
  expect_false(is.na(d$label[d$metric == "ai"]))
  expect_false(identical(d$label[d$metric == "ai"], "ai")) ## nrvtools' plain name, not the code
})

test_that("stand rows run biomass, age classes young to old, forested share, then leading types by label", {
  m <- c("lead_B", "age_Old", "forested", "lead_A", "age_Young", "biomassTg", "age_Mature", "age_Immature")
  cmp <- data.frame(layer = "L1", poly = "E", metric = m, group = "stand",
                    label = c("B leading", "Old", "Forested", "A leading", "Young", "Biomass", "Mature", "Immature"),
                    lo = 0, hi = 1, mid = 0.5, current = 0.5, position = "inside")
  cmp <- rbind(cmp, transform(cmpRows("E")[1:2, ], group = "landscape"))
  d <- ccFn(".currentConditionAllData", cmp)
  expect_identical(d$metric[d$group == "stand"], c("biomassTg", "age_Young", "age_Immature", "age_Mature", "age_Old",
                                                   "forested", "lead_A", "lead_B"))
  expect_identical(d$metric[d$group == "landscape"], c("ai", "area_mn")) ## unchanged
})

test_that("age classes are ordered young to old, in the stand overview too", {
  m <- c("lead_A", "age_Old", "forested", "age_Young", "biomassTg", "age_Mature", "age_Immature")
  cmp <- data.frame(poly = "E", metric = m, lo = 0, hi = 1, mid = 0.5, current = 0.5, position = "inside")
  pages <- ccFn(".standOverviewPages", cmp)
  expect_identical(
    as.character(unique(pages[[1]]$panel)),
    ccFn(".standLabel", c("biomassTg", "age_Young", "age_Immature", "age_Mature", "age_Old", "forested", "lead_A"))
  )
})

test_that("unit names 'A B' and 'A_B' give two files, whatever the order they come in", {
  nm <- function(x) ccFn(".ccAllNames", x)
  a <- nm(c("A B", "A_B"))
  expect_length(unique(tolower(a)), 2L)
  expect_true(all(startsWith(a, "A_B"))) ## still readable
  expect_identical(nm(c("A_B", "A B"))[c(2, 1)], a) ## the same file for each unit
  expect_length(unique(tolower(nm(c("ab", "AB", "a/b", "a b")))), 4L) ## also on a case-insensitive disk
  expect_length(unique(nm(c("A B", "A_B", "A_B-2"))), 3L) ## nor does a suffix collide with a unit
  expect_identical(unname(nm(c("ELF", "12.3", "a/b c"))), c("ELF", "12.3", "a_b_c")) ## others keep their name
  expect_identical(unname(nm("")), "unit")
})

## the sim and the module's `P` for the event: outside an event there is no current module
ccEvent <- function(comparisons, events = c("lm", "st", "cc"), .env = parent.frame()) {
  withr::local_package("SpaDES.core", .local_envir = .env)
  mod <- new.env()
  mod$nrvComparisons <- comparisons
  fn <- moduleFn("currentConditionAll", mod)
  environment(fn)$P <- function(sim) SpaDES.core::params(sim)$NRV_summary
  sim <- SpaDES.core::simInit(times = list(start = 0, end = 10),
                              params = list(NRV_summary = list(mode = "multi", reps = 1L, .plots = "png",
                                                               postprocessEvents = events)),
                              modules = moduleName, objects = list(reportingPolygons = list(ELF = "placeholder")),
                              paths = testPaths)
  list(sim = sim, fn = fn, dir = file.path(dirname(SpaDES.core::outputPath(sim)), "postprocess", "figures", "all"))
}

test_that("same-named units in two layers give two figures, registered as outputs", {
  e <- ccEvent(list(a = cmpRows("ELF", "L1"), b = cmpRows("ELF", "L2"), c = cmpRows("A B", "L2"),
                    d = cmpRows("A_B", "L2")))
  sim <- e$fn(e$sim)
  f <- file.path(e$dir, c(file.path("L1", "ELF.png"), file.path("L2", "ELF.png")))
  expect_true(all(file.exists(f)))
  expect_length(list.files(file.path(e$dir, "L2"), "\\.png$"), 3L) ## ELF, 'A B' and 'A_B'
  expect_length(list.files(file.path(e$dir, "L1"), "\\.png$"), 1L)
  expect_true(all(list.files(e$dir, "\\.png$", recursive = TRUE, full.names = TRUE) %in%
                    normalizePath(SpaDES.core::outputs(sim)$file)))
})

test_that("the cc event empties figures/all/ before it draws, and the other events leave it alone", {
  e <- ccEvent(list(a = cmpRows("ELF", "L1")))
  dir.create(file.path(e$dir, "Gone"), recursive = TRUE)
  stale <- file.path(e$dir, "Gone", "old.png")
  file.create(stale)
  e$fn(e$sim)
  expect_false(file.exists(stale))
  expect_true(file.exists(file.path(e$dir, "L1", "ELF.png")))
  clear <- moduleFn(".clearFigDirs", new.env())
  stale <- file.path(e$dir, "L1", "stale.png")
  file.create(stale)
  clear(e$sim, c("lm", "st"))
  expect_true(file.exists(stale))
  clear(e$sim, "cc")
  expect_false(file.exists(stale))
})
