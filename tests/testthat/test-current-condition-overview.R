## The current-condition overview: where each metric's year-0 value sits against its NRV range, on the
## metric's real scale. Ported from the PE fork (PE#12); drawn by the opt-in `stability` event.

## two metrics on very different scales, with a transient at time 100 before the later times
overviewEnv <- function(polys = "ELF") {
  g <- expand.grid(poly = polys, time = c(100, 200, 300), metric = c("ai", "area_mn"), stringsAsFactors = FALSE)
  base <- ifelse(g$metric == "ai", 10, 5000) * ifelse(g$time == 100, 5, 1)
  data.frame(poly = g$poly, time = g$time, metric = g$metric, mean = base, min = base - 2, max = base + 2)
}
overviewCC <- function(polys = "ELF") {
  g <- expand.grid(poly = polys, metric = c("ai", "area_mn"), stringsAsFactors = FALSE)
  data.frame(poly = g$poly, metric = g$metric, mean = ifelse(g$metric == "ai", 11, 5000))
}

test_that("the NRV range spans the NRV years, not the times before them", {
  nrvRange <- moduleFn(".nrvRange", new.env())
  r <- nrvRange(overviewEnv(), period = c(100, 300))
  r <- r[r$metric == "ai", ]
  expect_equal(c(r$lo, r$hi), c(8, 52)) ## the time-100 envelope is part of the NRV
  expect_equal(r$mid, median(c(50, 10, 10)))
  r2 <- nrvRange(overviewEnv(), period = c(200, 300))
  expect_equal(c(r2$lo[r2$metric == "ai"], r2$hi[r2$metric == "ai"]), c(8, 12))
})

test_that("an all-NA year does not blank a subregion's NRV range", {
  nrvRange <- moduleFn(".nrvRange", new.env())
  env <- overviewEnv()
  i <- env$metric == "ai" & env$time == 200
  env[i, c("mean", "min", "max")] <- NA_real_
  r <- nrvRange(env, period = c(100, 300))
  r <- r[r$metric == "ai", ]
  expect_equal(c(r$lo, r$hi), c(8, 52))
  expect_equal(r$mid, median(c(50, 10)))
  ## a metric with no value at all in the period has no range, and no warning
  env$min[env$metric == "area_mn"] <- NA_real_
  env$max[env$metric == "area_mn"] <- NA_real_
  env$mean[env$metric == "area_mn"] <- NA_real_
  expect_no_warning(r <- nrvRange(env, period = c(100, 300)))
  expect_false("area_mn" %in% r$metric)
})

test_that("the current condition is classed inside, below or above the NRV range", {
  ccData <- moduleFn(".ccOverviewData", new.env())
  cc <- data.frame(poly = "ELF", metric = c("ai", "area_mn"), mean = c(2, 90000))
  d <- ccData(overviewEnv(), cc, c(100, 300))
  expect_equal(d$position[d$metric == "ai"], "below")
  expect_equal(d$position[d$metric == "area_mn"], "above")
  cc$mean <- c(11, 5000)
  expect_setequal(ccData(overviewEnv(), cc, c(100, 300))$position, "inside")
})

test_that("panels are labelled with nrvtools' metric titles, with the subregion", {
  ccData <- moduleFn(".ccOverviewData", new.env())
  d <- ccData(overviewEnv(c("A", "B")), overviewCC(c("A", "B")), c(100, 300))
  title <- moduleFn(".metricTitle", new.env())
  mets <- c("ai", "area_mn")
  want <- vapply(mets, title, character(1), level = "landscape")
  expect_setequal(d$panel, paste0(rep(c("A", "B"), each = 2), ": ", rep(want, 2)))
})

test_that("a 20-unit layer gives pages of at most 24 panels", {
  ccData <- moduleFn(".ccOverviewData", new.env())
  pages <- moduleFn(".ccOverviewPages", new.env())
  polys <- sprintf("P%02d", 1:20)
  d <- ccData(overviewEnv(polys), overviewCC(polys), c(100, 300)) ## 20 x 2 = 40 panels
  pp <- pages(d)
  expect_equal(length(pp), 2L)
  expect_lte(max(vapply(pp, function(p) length(unique(p$panel)), integer(1))), 24L)
  expect_equal(sort(unlist(lapply(pp, function(p) unique(p$panel)))), sort(unique(d$panel))) ## each panel once
  expect_equal(length(pages(d[d$poly == "P01", ])), 1L)
})

test_that("a page is drawn with only its own panels, on each metric's real scale", {
  ccData <- moduleFn(".ccOverviewData", new.env())
  pages <- moduleFn(".ccOverviewPages", new.env())
  draw <- moduleFn(".currentConditionOverview", new.env())
  polys <- sprintf("P%02d", 1:20)
  d <- ccData(overviewEnv(polys), overviewCC(polys), c(100, 300))
  pg <- pages(d)[[2L]]
  gg <- draw(pg, title = "t")
  b <- ggplot2::ggplot_build(gg)
  expect_s3_class(gg$facet, "FacetWrap")
  expect_true(gg$facet$params$free$x)
  expect_equal(length(unique(b$layout$layout$PANEL)), length(unique(pg$panel)))
  x <- unlist(lapply(b$data, function(l) c(l$x, l$xmin, l$xend)))
  expect_true(all(c(8, 52, 11, 5000) %in% x)) ## raw min, max and current condition; not 0-1
  expect_gt(max(x), 1000)
  src <- file.path(moduleRoot, paste0(moduleName, ".R"))
  skip_if_not(file.exists(src), "module source not found")
  expect_false("facet_wrap_paginate" %in% all.names(parse(src, keep.source = FALSE)))
})

test_that("the canvas follows the page's panel count, up to 24 panels", {
  lay <- moduleFn(".ccOverviewLayout", new.env())
  expect_equal(lay(24L)$nrow, 6L)
  expect_gt(lay(24L)$height, lay(12L)$height)
  expect_gt(lay(12L)$height, lay(4L)$height)
  expect_equal(lay(1L)$ncol, 1L)
})
