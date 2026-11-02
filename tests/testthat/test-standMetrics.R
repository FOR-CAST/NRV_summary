## The stand-level measures of the `st` event, from synthetic snapshot values.

test_that("a snapshot's stand-level row: age classes of forest, forested and leading shares of the area, biomass", {
  ## 6 reporting-area pixels: 4 forested (ages 10, 50, 90, 150), one with no age, one with pixel group 0
  sam <- c(10, 50, 90, 150, NA, 30)
  pgm <- c(1L, 1L, 2L, 2L, 3L, 0L)
  vtm <- c("Pice_mar", "Pice_mar", "Mixed", "Pice_mar", NA, NA)
  cd <- data.frame(pixelGroup = c(1L, 1L, 2L), B = c(100, 200, 1000)) ## g/m2: 300 in group 1, 1000 in group 2
  r <- .standRow(sam, vtm, pgm, cd, pixelHa = 6.25, ageClassCutOffs = c(0, 40, 80, 120),
                 ageClasses = c("Young", "Immature", "Mature", "Old"))
  v <- setNames(r$value, r$metric)
  expect_equal(v[paste0("age_", c("Young", "Immature", "Mature", "Old"))], setNames(rep(25, 4), paste0("age_", c("Young", "Immature", "Mature", "Old"))))
  expect_equal(v[["forested"]], 100 * 4 / 6)
  expect_equal(v[["lead_Pice_mar"]], 100 * 3 / 6)
  expect_equal(v[["lead_Mixed"]], 100 * 1 / 6)
  ## B over the pixels: 2 x 300 + 2 x 1000 g/m2, x 6.25 ha x 1e4 m2/ha, in Tg
  expect_equal(v[["biomassTg"]], (2 * 300 + 2 * 1000) * 6.25 * 1e4 / 1e12)
})

test_that("snapshot files pair by rep folder and year", {
  f <- function(kind, ext) c(sprintf("o/rep2/%s_year2120.%s", kind, ext), sprintf("o/rep1/%s_year2020.%s", kind, ext),
                             sprintf("o/rep1/%s_year2120.%s", kind, ext))
  s <- .snapshotTable(f("standAgeMap", "tif"), f("vegTypeMap", "tif"), f("pixelGroupMap", "tif"), f("cohortData", "qs2"))
  expect_identical(nrow(s), 3L)
  i <- s$rep == "rep2"
  expect_identical(s$time[i], 2120)
  expect_identical(s$cd[i], "o/rep2/cohortData_year2120.qs2")
  expect_identical(s$vtm[i], "o/rep2/vegTypeMap_year2120.tif")
})

test_that("a type absent from a snapshot leads 0% of it, also one found only in the current condition", {
  rows <- data.frame(rep = c(1, 1, 2), time = c(0, 100, 0), metric = c("lead_A", "lead_B", "lead_A"), value = c(5, 7, 3))
  out <- .fillAbsentMetrics(rows, "lead_C")
  expect_identical(nrow(out), 2L * 2L * 3L)
  expect_equal(out$value[out$rep == 2 & out$time == 100 & out$metric == "lead_A"], 0)
  expect_equal(sum(out$value[out$metric == "lead_C"]), 0)
})

test_that("the current condition is inside, above or below the NRV range of the NRV years only", {
  env <- data.frame(poly = "E", metric = rep(c("a", "b", "c"), each = 3), time = rep(c(0, 700, 1000), 3),
                    mean = 5, min = c(0, 2, 3, 0, 2, 3, 0, 2, 3), max = c(100, 8, 9, 100, 8, 9, 100, 8, 9))
  cc <- data.frame(poly = "E", metric = c("a", "b", "c"), mean = c(9, 50, 1))
  r <- .nrvComparison(env, cc, c(700, 1000))
  r <- r[order(r$metric), ]
  expect_identical(r$position, c("inside", "above", "below")) ## 9 is the top of the range: inside
  expect_equal(r$lo, rep(2, 3))
  expect_equal(r$hi, rep(9, 3))
})

test_that("stand-level labels name the leading type by sppEquiv", {
  se <- data.frame(LandR = "Pice_mar", EN_generic_full = "Black spruce")
  expect_identical(.standLabel(c("lead_Pice_mar", "age_Old", "forested", "biomassTg"), se),
                   c("Black spruce leading, share of the area (%)", "Old forest, share of forested area (%)",
                     "Forested share of the area (%)", "Total tree biomass (Tg)"))
})

test_that("the stand-level figure draws reps, their median and the current condition", {
  s <- data.frame(rep = rep(1:2, each = 3), time = rep(c(0, 500, 1000), 2), label = "Old", value = 1:6)
  g <- .standSeriesPanels(list(series = s, cc = data.frame(time = 0, label = "Old", value = 2),
                               period = c(700, 1000), title = "t", ylab = "y"))
  expect_s3_class(g, "ggplot")
  expect_silent(ggplot2::ggplot_build(g))
})
