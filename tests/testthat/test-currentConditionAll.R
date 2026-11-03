## One figure with every NRV metric of every group, each metric's NRV range rescaled to 0-1.

cmpRows <- function(poly = "ELF") {
  data.frame(
    poly = poly, metric = c("ai", "area_mn", "age_Old", "lead_A", "lead_Zero", "forested"),
    group = c("landscape", "landscape", "stand", "stand", "stand", "stand"),
    label = c("Aggregation index", "Mean patch area", "Old forest", "A leading", "Zero leading", "Forested"),
    lo = c(10, 100, 40, 1, 0, 50), hi = c(20, 200, 60, 3, 0, 50), mid = c(15, 150, 50, 2, 0, 50),
    current = c(15, 50, 90, 2, 0, 50), position = c("inside", "below", "above", "inside", "inside", "inside")
  )
}

test_that("a value is rescaled so the NRV range is 0 to 1", {
  expect_equal(.nrvScale(c(10, 15, 20, 5, 30), 10, 20), c(0, 0.5, 1, -0.5, 2))
})

test_that("a range of one value gets a fixed small width around it", {
  x <- .nrvScale(c(49.5, 50, 50.5), 50, 50)
  expect_equal(x[2], 0.5) ## the value is the middle of the bar
  expect_equal(.nrvScale(c(0, 1), 0, 0), c(0.5, 100.5))
  expect_true(all(is.finite(x)))
})

test_that("rows of every group appear, and a leading type that is 0 everywhere is dropped", {
  d <- .currentConditionAllData(cmpRows())
  expect_setequal(d$metric, c("ai", "area_mn", "age_Old", "lead_A", "forested"))
  expect_setequal(d$group, c("landscape", "stand"))
  expect_false("lead_Zero" %in% d$metric)
  expect_equal(d$xCurrent[d$metric == "age_Old"], 2.5) ## (90 - 40) / 20 ...
})

test_that("far-out points are clipped to the axis limit and keep their true value", {
  d <- .currentConditionAllData(cmpRows())
  expect_true(all(d$xPlot >= .ccAllLimits[1] & d$xPlot <= .ccAllLimits[2]))
  expect_equal(d$xPlot[d$metric == "age_Old"], .ccAllLimits[2])
  expect_true(d$clipped[d$metric == "age_Old"])
  expect_false(d$clipped[d$metric == "ai"])
})

test_that("the figure is one panel, with the metrics as rows coloured by position", {
  gg <- .currentConditionAll(.currentConditionAllData(cmpRows()))
  b <- ggplot2::ggplot_build(gg)
  expect_s3_class(gg$facet, "FacetNull")
  expect_equal(length(unique(b$layout$layout$PANEL)), 1L)
  expect_identical(gg$labels$x, .ccAllXlab)
  expect_true(all(unname(.positionColour) %in% unlist(lapply(b$data, function(l) l$colour)))) ## inside, below, above
})

test_that("the file name carries the polygon only when there are several", {
  expect_identical(.ccAllName("ELF", 1L), "current_condition_all")
  expect_identical(.ccAllName("12.3", 2L), "current_condition_all_12.3")
  expect_identical(.ccAllName("a/b c", 2L), "current_condition_all_a_b_c")
})

test_that("each group's comparison is recorded, with its label, for the combined figure", {
  mod <- new.env()
  record <- moduleFn(".recordComparison", mod)
  record("lm", data.frame(poly = "ELF", metric = "ai", lo = 1, hi = 2, mid = 1.5, current = 3, position = "above"))
  record("st", data.frame(poly = "ELF", metric = "age_Old", lo = 1, hi = 2, mid = 1.5, current = 1, position = "inside",
                          label = "Old forest"))
  d <- do.call(rbind, mod$nrvComparisons)
  expect_setequal(d$group, c("landscape", "stand"))
  expect_equal(d$label[d$metric == "ai"], "Aggregation index")
  expect_equal(d$label[d$metric == "age_Old"], "Old forest")
})

test_that("the combined figure is written once per polygon and registered as an output", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  mod$nrvComparisons <- list(a = cmpRows("ELF"))
  fn <- moduleFn("currentConditionAll", mod)
  environment(fn)$P <- function(sim) SpaDES.core::params(sim)$NRV_summary ## outside an event there is no current module
  sim <- SpaDES.core::simInit(times = list(start = 0, end = 10),
                              params = list(NRV_summary = list(mode = "multi", reps = 1L, .plots = "png")),
                              modules = moduleName, objects = list(reportingPolygons = list(ELF = "placeholder")),
                              paths = testPaths)
  sim <- fn(sim)
  f <- file.path(dirname(SpaDES.core::outputPath(sim)), "postprocess", "figures", "all", "current_condition_all.png")
  expect_true(file.exists(f))
  expect_true(normalizePath(f) %in% normalizePath(SpaDES.core::outputs(sim)$file))
})

test_that("the NRV mid is the median over the NRV years, not the mean", {
  env <- data.frame(poly = "E", metric = "a", time = c(0, 700, 800, 1000), mean = c(1000, 1, 2, 9), min = 0, max = 10)
  expect_equal(.nrvRange(env, c(700, 1000))$mid, 2) ## the mean would be 4; time 0 is outside the NRV years
})

test_that("the NRV tick is labelled median on the combined and the per-group figures", {
  expect_identical(.nrvMidLabel, "NRV median")
  cc <- data.frame(poly = "ELF", metric = "ai", mean = 11)
  env <- data.frame(time = c(100, 200), poly = "ELF", metric = "ai", mean = 10, min = 8, max = 12)
  g1 <- ggplot2::ggplot_build(.currentConditionOverview(list(env = env, cc = cc, period = c(100, 200))))
  g2 <- ggplot2::ggplot_build(.currentConditionAll(.currentConditionAllData(cmpRows())))
  expect_identical(g1$plot$scales$get_scales("shape")$get_labels(), "NRV median")
  expect_identical(g2$plot$scales$get_scales("shape")$get_labels(), "NRV median")
})

test_that("stand rows run biomass, age classes young to old, forested share, then leading types by label", {
  m <- c("lead_B", "age_Old", "forested", "lead_A", "age_Young", "biomassTg", "age_Mature", "age_Immature")
  cmp <- data.frame(poly = "E", metric = m, group = "stand", label = c("B leading", "Old", "Forested", "A leading", "Young",
                    "Biomass", "Mature", "Immature"), lo = 0, hi = 1, mid = 0.5, current = 0.5, position = "inside")
  cmp <- rbind(cmp, transform(cmpRows("E")[1:2, ], group = "landscape"))
  d <- .currentConditionAllData(cmp)
  expect_identical(d$metric[d$group == "stand"], c("biomassTg", "age_Young", "age_Immature", "age_Mature", "age_Old",
                                                   "forested", "lead_A", "lead_B"))
  expect_identical(d$metric[d$group == "landscape"], c("ai", "area_mn")) ## unchanged
})
