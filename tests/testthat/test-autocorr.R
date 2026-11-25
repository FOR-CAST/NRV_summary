## The autocorrelation time of the landscape gives the interval at which summary snapshots are independent.
## Ported from the PE fork (PE#14, PE#16); the checks are opt-in here (`recordAnnualSeries`, the `autocorr` event).

## the module functions under test, looked up when called (see moduleFn())
.acf <- function(...) moduleFn(".acf", new.env())(...)
.autocorrTime <- function(...) moduleFn(".autocorrTime", new.env())(...)
.firstLagBelow <- function(...) moduleFn(".firstLagBelow", new.env())(...)
.essAR1 <- function(...) moduleFn(".essAR1", new.env())(...)
.autocorrSummary <- function(...) moduleFn(".autocorrSummary", new.env())(...)
.annualSeriesRow <- function(...) moduleFn(".annualSeriesRow", new.env())(...)
.bindAnnualRows <- function(...) moduleFn(".bindAnnualRows", new.env())(...)
.annualSeriesLabel <- function(...) moduleFn(".annualSeriesLabel", new.env())(...)

ar1 <- function(phi, n = 20000, seed = 1) {
  set.seed(seed)
  as.numeric(stats::arima.sim(list(ar = phi), n = n))
}

test_that("tau of an AR(1) series with phi = 0.9 is (1 + phi) / (1 - phi) = 19", {
  tau <- .autocorrTime(.acf(ar1(0.9, n = 2e5), 200))
  expect_equal(tau, 19, tolerance = 0.1)
})

test_that("white noise has tau of about 1", {
  set.seed(2)
  expect_equal(.autocorrTime(.acf(stats::rnorm(2e5), 200)), 1, tolerance = 0.1)
})

test_that("the first lag with ACF below 0.1 is where it is, NA when it never is", {
  expect_identical(.firstLagBelow(c(1, 0.5, 0.2, 0.05, 0)), 3L)
  expect_identical(.firstLagBelow(c(1, 0.5, 0.4)), NA_integer_)
})

## reps x years of AR(1) (or other) series, as the summary reads them
repTables <- function(f, reps = 5, years = 300, seed = 1) {
  lapply(seq_len(reps), function(r) {
    set.seed(seed * 100 + r)
    data.frame(time = seq_len(years), x = f(years))
  })
}
nrv <- c(1, 300)

## 15 reps x 301 NRV years of AR(1) with phi 0.95 (tau 39), as LandWeb's NRV: with each rep centred on its own
## mean, the estimated tau is about 0.75 of the truth (the sample ACF is biased toward 0 when tau is a sizeable
## fraction of the series); around the replicates' common mean it is about 0.95. 30 fixed seeds averaged.
test_that("15 reps x 301 years of AR(1) with phi 0.95 give a mean tau within 15% of 39", {
  tau <- vapply(1:30, function(seed) {
    .autocorrSummary(repTables(function(n) as.numeric(stats::arima.sim(list(ar = 0.95), n)), reps = 15,
                               years = 301, seed = seed), c(1, 301), 200, 50)$fit$tau
  }, numeric(1))
  expect_gt(mean(tau) / 39, 0.85)
  expect_lt(mean(tau) / 39, 1.15)
})

test_that(".acf centres on the series' own mean unless given a centre and a scale", {
  set.seed(3)
  x <- as.numeric(stats::arima.sim(list(ar = 0.7), 200))
  expect_equal(.acf(x, 20), as.numeric(stats::acf(x, lag.max = 20, plot = FALSE)$acf))
  ## around a given centre, with a given sum of squares: the lagged products over that sum
  d <- x - 1
  expect_equal(.acf(x, 5, centre = 1, scale = 2 * sum(d^2))[3], sum(d[1:198] * d[3:200]) / (2 * sum(d^2)))
})

test_that("the ACF is pooled around the replicates' common mean", {
  ## two reps with different levels and no variation inside: around their own means there is nothing to
  ## correlate, around the common mean each is perfectly persistent
  tables <- list(data.frame(time = 1:60, x = 0), data.frame(time = 1:60, x = 1))
  fit <- .autocorrSummary(tables, c(1, 60), 20, 10)$fit
  expect_gt(fit$phi, 0.9)
})

test_that("each replicate's ACF is a correlation (1 at lag 0, within [-1, 1]); the pooled one is the pooled estimator", {
  ## replicates of very different variance around one level: the per-rep curves are scaled by their own sum of squares
  set.seed(11)
  sd <- c(1, 5, 0.2)
  tables <- lapply(sd, function(k) data.frame(time = 1:200, x = k * as.numeric(stats::arima.sim(list(ar = 0.8), 200))))
  res <- .autocorrSummary(tables, c(1, 200), 30, 10)
  perRep <- res$acf[res$acf$rep != "pooled", ]
  expect_equal(perRep$acf[perRep$lag == 0], rep(1, 3))
  expect_true(all(abs(perRep$acf) <= 1 + 1e-12))
  ## the pooled curve: summed lagged products over the summed squares, around the common mean
  x <- lapply(tables, `[[`, "x")
  m <- mean(unlist(x))
  pooledAt <- function(k) {
    sum(vapply(x, function(v) sum((v[seq_len(200 - k)] - m) * (v[seq_len(200 - k) + k] - m)), 0)) /
      sum(vapply(x, function(v) sum((v - m)^2), 0))
  }
  pooled <- res$acf[res$acf$rep == "pooled", ]
  expect_equal(pooled$acf[pooled$lag %in% 0:5], vapply(0:5, pooledAt, 0))
  expect_lte(res$fit$rMax, 1)
  expect_gte(res$fit$rMin, -1)
})

test_that("a series with fewer than 3 values in a replicate is 'tooFew', not 'constant'", {
  tables <- list(data.frame(time = 1:100, a = stats::rnorm(100), b = 1, c = stats::rnorm(100)),
                 data.frame(time = 1:100, a = c(1, 2, rep(NA, 98)), b = 1, c = stats::rnorm(100)))
  res <- .autocorrSummary(tables, c(1, 100), 20, 10)
  expect_identical(res$tooFew, "a")
  expect_identical(res$constant, "b")
  expect_setequal(res$undefined, c("a", "b"))
  expect_false(is.na(res$limiting))
})

test_that("a merged species code is labelled with all of its names", {
  sppEquiv <- data.frame(LandR = c("Pice_mar", "Pice_gla", "Pice_gla", "Betu_pap"),
                         EN_generic_full = c("Black spruce", "White spruce", "Engelmann spruce", "Paper birch"))
  expect_identical(.annualSeriesLabel(c("lead_Pice_gla", "lead_Betu_pap"), sppEquiv),
                   c("Leading: White spruce / Engelmann spruce", "Leading: Paper birch"))
  ## the same name twice is said once
  dup <- data.frame(LandR = c("Pice_gla", "Pice_gla"), EN_generic_full = c("Spruce", "Spruce"))
  expect_identical(.annualSeriesLabel("lead_Pice_gla", dup), "Leading: Spruce")
})

test_that("r at the summary interval is close to phi^interval and its reps bracket it", {
  tables <- repTables(function(n) as.numeric(stats::arima.sim(list(ar = 0.95), n)), reps = 20, years = 2000)
  at100 <- .autocorrSummary(tables, c(1, 2000), 200, 100)$fit
  expect_lt(abs(at100$r - 0.95^100), 0.05)
  at20 <- .autocorrSummary(tables, c(1, 2000), 200, 20)$fit
  expect_equal(at20$r, 0.95^20, tolerance = 0.1)
  expect_lte(at20$rMin, at20$r)
  expect_gte(at20$rMax, at20$r)
})

test_that("the effective snapshots follow n (1 - r) / (1 + r), with r clipped to [0, 1)", {
  expect_equal(.essAR1(10, 0.5), 10 / 3)
  expect_equal(.essAR1(10, -0.3), 10)
  expect_lt(.essAR1(10, 1.2), 0.1 * 10)
})

test_that("the misfit is near 0 for AR(1) and large for an oscillating AR(2) series", {
  ar1Fit <- .autocorrSummary(repTables(function(n) as.numeric(stats::arima.sim(list(ar = 0.9), n)), years = 2000),
                             c(1, 2000), 100, 100)$fit
  osc <- .autocorrSummary(repTables(function(n) as.numeric(stats::arima.sim(list(ar = c(1.4, -0.8)), n)), years = 2000),
                          c(1, 2000), 100, 100)$fit
  expect_lt(ar1Fit$ar1Misfit, 0.1)
  expect_gt(osc$ar1Misfit, 0.3)
})

test_that("the recommended thin is the ceiling of the largest tau; the Geyer tau is kept; a constant series is NA", {
  n <- 600
  tables <- lapply(1:3, function(r) {
    set.seed(r)
    data.frame(time = seq_len(n), slow = ar1(0.9, n, r), fast = ar1(0.3, n, 10 + r), constant = 1)
  })
  res <- .autocorrSummary(tables, c(1, n), maxLag = 100, interval = 100)
  fit <- res$fit
  expect_identical(res$thin, ceiling(max(fit$tau, na.rm = TRUE)))
  expect_identical(res$limiting, "slow")
  expect_true(is.na(fit$tau[fit$series == "constant"]))
  expect_true(all(c("phi", "tau", "ar1Misfit", "r", "rMin", "rMax", "effectiveSnapshots", "tauGeyer") %in% names(fit)))
  expect_equal(res$nSnapshots, 3 * (floor((n - 1) / 100) + 1))
  expect_setequal(unique(res$acf$rep), c("1", "2", "3", "pooled"))
})

test_that("the annual recorder returns one row per call, equal to the same scalars from the maps", {
  ## 4 x 4 pixels of 100 m; pixel group 3 has no cohorts (not forest)
  pgm <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 400, ymin = 0, ymax = 400, vals = rep(c(1, 2, 3, 1), 4))
  names(pgm) <- "pixelGroup"
  cd <- data.table::data.table(pixelGroup = c(1L, 1L, 2L, 2L), speciesCode = c("Pice_mar", "Betu_pap", "Betu_pap", "Pice_mar"),
                               age = c(30L, 60L, 150L, 160L), B = c(1000L, 500L, 2000L, 100L))
  sppEquiv <- data.table::data.table(LandR = c("Pice_mar", "Betu_pap"), Type = c("Conifer", "Deciduous"),
                                     EN_generic_full = c("Black spruce", "Paper birch"))
  cuts <- c(0L, 40L, 80L, 120L)
  row <- function(t) .annualSeriesRow(t, cd, tabulate(terra::values(pgm, mat = FALSE)), 1, cuts, 0.8, 2L, sppEquiv, "LandR")
  rows <- lapply(c(5, 6), row)
  expect_identical(nrow(rows[[1]]), 1L)
  expect_identical(rows[[1]]$time, 5)
  expect_identical(nrow(.bindAnnualRows(rows)), 2L)

  ## the same scalars from the maps map_generators makes
  sam <- LandR::standAgeMapGenerator(cd, pgm, weight = "biomass", doAssertion = FALSE)
  vtm <- LandR::vegTypeMapGenerator(cd, pgm, 0.8, mixedType = 2L, sppEquiv = sppEquiv, sppEquivCol = "LandR",
                                    colors = c(Pice_mar = "black", Betu_pap = "white", Mixed = "grey"), doAssertion = FALSE)
  age <- terra::values(sam, mat = FALSE)
  forest <- !is.na(age)
  lv <- terra::levels(vtm)[[1L]]
  lead <- tabulate(match(terra::values(vtm, mat = FALSE)[forest], lv[[1L]]), nbins = nrow(lv)) / sum(forest)
  expected <- c(propYoung = mean(age[forest] < 40), propOld = mean(age[forest] >= 120), meanStandAge = mean(age[forest]),
                setNames(lead, paste0("lead_", lv[[2L]])))
  got <- unlist(rows[[1]][setdiff(names(expected), "totalBiomassTg")])
  expect_equal(got[order(names(got))], expected[order(names(expected))])
  ## biomass: 8 pixels of group 1 (1500 g/m2) and 4 of group 2 (2100 g/m2), 1 ha each
  expect_equal(rows[[1]]$totalBiomassTg, (8 * 1500 + 4 * 2100) * 0.01 / 1e6)
})

test_that("the annual row is the same as with a join of the leading species onto every cohort", {
  ## pixelGroups 1 and 2 have two deciduous species (summed for the mixedwood decision), 3 is pure conifer,
  ## 4 has a tie between two conifers
  cd <- data.table::data.table(
    pixelGroup = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 4L, 4L),
    speciesCode = factor(c("Pice_mar", "Betu_pap", "Popu_tre", "Pice_mar", "Betu_pap", "Popu_tre", "Pice_mar",
                           "Pinu_ban", "Pice_mar")),
    age = c(30L, 60L, 90L, 50L, 150L, 160L, 70L, 20L, 80L),
    B = c(2500L, 1500L, 1000L, 7000L, 2000L, 1000L, 800L, 500L, 500L)
  )
  sppEquiv <- data.table::data.table(LandR = c("Pice_mar", "Betu_pap", "Popu_tre", "Pinu_ban"),
                                     Type = c("Conifer", "Deciduous", "Deciduous", "Conifer"))
  pgm <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 400, ymin = 0, ymax = 400, vals = rep(c(1, 2, 3, 4), 4))
  nPerGroup <- tabulate(terra::values(pgm, mat = FALSE))
  cuts <- c(0L, 40L, 80L, 120L)
  got <- .annualSeriesRow(7, cd, nPerGroup, 1, cuts, 0.75, 2L, sppEquiv, "LandR")

  ## the row as it was made: vegTypeGenerator() joins `leading` onto every cohort row, one row per group kept
  lead <- LandR::vegTypeGenerator(cd, 0.75, mixedType = 2L, sppEquiv = sppEquiv, sppEquivCol = "LandR",
                                  doAssertion = FALSE)
  lead <- lead[!duplicated(lead$pixelGroup), ]
  expect_identical(sort(as.character(lead$leading[order(lead$pixelGroup)])),
                   sort(c("Mixed", "Mixed", "Pice_mar", "Pice_mar")))  ## (the case is not trivial)
  byLead <- tapply(nPerGroup[lead$pixelGroup], as.character(lead$leading), sum) / sum(nPerGroup[lead$pixelGroup])
  expect_identical(unlist(got[paste0("lead_", names(byLead))], use.names = FALSE), as.numeric(byLead))
  expect_setequal(grep("^lead_", names(got), value = TRUE), paste0("lead_", names(byLead)))
  expect_identical(.annualSeriesRow(7, cd, nPerGroup, 1, cuts, 0.75, 2L, sppEquiv, "LandR"), got)
})

test_that("series have full names, a species by its full name, one label per series", {
  sppEquiv <- data.frame(LandR = c("Pice_mar", "Betu_pap"), EN_generic_full = c("Black spruce", "Paper birch"))
  series <- c("propYoung", "lead_Pice_mar", "lead_Betu_pap", "lead_Mixed", "meanStandAge")
  expect_identical(.annualSeriesLabel(series, sppEquiv),
                   c("Proportion of forest young", "Leading: Black spruce", "Leading: Paper birch",
                     "Leading: Mixed", "Mean stand age"))
  expect_identical(.annualSeriesLabel(series[1:2]), c("Proportion of forest young", "Leading: Pice_mar"))
})

test_that("lead_ labels fall back to the species code without EN_generic_full", {
  noFull <- data.frame(LandR = c("Pice_mar", "Betu_pap"), Type = c("Conifer", "Deciduous"))
  series <- c("propYoung", "lead_Pice_mar", "lead_Mixed")
  expect_identical(.annualSeriesLabel(series, noFull),
                   c("Proportion of forest young", "Leading: Pice_mar", "Leading: Mixed"))
})

test_that("tau is finite when the series is shorter than the maximum lag", {
  ## the NRV years of nrvWindow = 0.3 can be fewer than autocorrMaxLag (200)
  tables <- repTables(function(n) as.numeric(stats::arima.sim(list(ar = 0.8), n)), reps = 4, years = 60)
  res <- .autocorrSummary(tables, c(1, 60), maxLag = 200, interval = 50)
  fit <- res$fit
  expect_true(all(is.finite(c(fit$phi, fit$tau, fit$tauGeyer, fit$ar1Misfit, fit$r, fit$effectiveSnapshots))))
  expect_true(is.finite(res$thin))
  ## no lag beyond the series' length is made up
  expect_lte(max(res$acf$lag), 59L)
  ## an interval longer than the series has no correlation to report
  long <- .autocorrSummary(tables, c(1, 60), maxLag = 200, interval = 100)$fit
  expect_true(is.na(long$r))
  expect_true(is.finite(long$tau))
})
