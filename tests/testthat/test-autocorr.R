## The annual landscape series that the `autocorr` event reads: the recorder and the series' labels. Ported from the
## PE fork (PE#14, PE#16); the checks are opt-in here (`recordAnnualSeries`, the `autocorr` event). The autocorrelation
## itself is nrvtools::autocorr_summary()'s, tested in nrvtools; the event is tested in test-annual-series.R.

## the module functions under test, looked up when called (see moduleFn())
.annualSeriesRow <- function(...) moduleFn(".annualSeriesRow", new.env())(...)
.bindAnnualRows <- function(...) moduleFn(".bindAnnualRows", new.env())(...)
.annualSeriesLabel <- function(...) moduleFn(".annualSeriesLabel", new.env())(...)

test_that("a merged species code is labelled with all of its names", {
  sppEquiv <- data.frame(LandR = c("Pice_mar", "Pice_gla", "Pice_gla", "Betu_pap"),
                         EN_generic_full = c("Black spruce", "White spruce", "Engelmann spruce", "Paper birch"))
  expect_identical(.annualSeriesLabel(c("lead_Pice_gla", "lead_Betu_pap"), sppEquiv),
                   c("Leading: White spruce / Engelmann spruce", "Leading: Paper birch"))
  ## the same name twice is said once
  dup <- data.frame(LandR = c("Pice_gla", "Pice_gla"), EN_generic_full = c("Spruce", "Spruce"))
  expect_identical(.annualSeriesLabel("lead_Pice_gla", dup), "Leading: Spruce")
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
