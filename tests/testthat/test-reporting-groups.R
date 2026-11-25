## Reporting groups (`sppEquivColReporting`): mode = "multi" analyses vegetation-type maps of groups of
## the simulated species, rebuilt from the saved cohorts. LandR::vegTypeMapGenerator() sums biomass
## within a species code only, so a group's leading type needs the codes recoded first.

## the reporting table of a small landscape: two spruce-group species, two pines, two broadleaves
reportingEquiv <- function() {
  data.table::data.table(
    LandWeb = c("Pice_mar", "Lari_lar", "Pinu_ban", "Pinu_con", "Popu_tre", "Betu_pap", "Pice_gla"),
    LandWebReport = c("Bl_Spruce", "Bl_Spruce", "Pine", "Pine", "Decid", "Decid", "Wh_Spruce"),
    Type = c("Conifer", "Conifer", "Conifer", "Conifer", "Deciduous", "Deciduous", "Conifer")
  )
}

## four pixel groups, one pixel each:
##  1: black spruce 40% + tamarack 40% + aspen 20%      -> Bl_Spruce (80% as a group)
##  2: aspen 50% + birch 30% + white spruce 20%          -> Decid (broadleaf 80%, not mixed)
##  3: jack pine 30% + lodgepole pine 30% + black spruce 40% -> Pine by group, black spruce by species
##  4: aspen 50% + black spruce 50%                      -> Mixed
reportingCohorts <- function() {
  data.table::data.table(
    pixelGroup = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L),
    speciesCode = c(
      "Pice_mar", "Lari_lar", "Popu_tre",
      "Popu_tre", "Betu_pap", "Pice_gla",
      "Pinu_ban", "Pinu_con", "Pice_mar",
      "Popu_tre", "Pice_mar"
    ),
    age = 50L,
    B = c(400L, 400L, 200L, 500L, 300L, 200L, 300L, 300L, 400L, 500L, 500L)
  )
}

reportingPGM <- function() {
  terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:3857",
              vals = 1:4)
}

## class label of each cell of a categorical SpatRaster
cellLabels <- function(r) {
  lv <- terra::levels(r)[[1L]]
  as.character(lv[[2L]][match(terra::values(r, mat = FALSE), lv[[1L]])])
}

test_that("the group map sums biomass within a group before deciding the leading type", {
  withr::local_package("data.table")
  grp <- LandWebUtils::sppEquiv_groups(reportingEquiv(), "LandWeb", "LandWebReport")
  colors <- moduleFn(".reportingColors", new.env())(NULL, grp)
  reportingVegTypeMap <- moduleFn(".reportingVegTypeMap", new.env())

  cd <- reportingCohorts()
  before <- data.table::copy(cd)
  vtm <- reportingVegTypeMap(cd, reportingPGM(), grp, colors, vegLeadingProportion = 0.75, mixedType = 2L)

  expect_identical(cellLabels(vtm), c("Bl_Spruce", "Decid", "Pine", "Mixed"))
  expect_identical(cd, before) ## the cohorts handed in are not recoded in place

  ## the species-level map, for contrast: black spruce leads pixel 3, no pine species does alone
  species <- LandR::vegTypeMapGenerator(
    reportingCohorts(), reportingPGM(), 0.75, mixedType = 2L,
    sppEquiv = reportingEquiv(), sppEquivCol = "LandWeb",
    colors = c(stats::setNames(rep("#000000", 7L), reportingEquiv()[["LandWeb"]]), Mixed = "#FFFFFF")
  )
  expect_identical(cellLabels(species)[[3L]], "Pice_mar")
})

test_that("colours default from LandR and must cover every group plus Mixed", {
  withr::local_package("data.table")
  grp <- LandWebUtils::sppEquiv_groups(reportingEquiv(), "LandWeb", "LandWebReport")
  reportingColors <- moduleFn(".reportingColors", new.env())

  expect_identical(names(reportingColors(NULL, grp)), c(grp$groups[["LandWebReport"]], "Mixed"))
  given <- c(Pine = "#E0A100", Decid = "#9BBF3B", Bl_Spruce = "#2A6FB0", Wh_Spruce = "#1B8A5A",
             Fir = "#19A7A0", Mixed = "#7D5BA6")
  expect_identical(names(reportingColors(given, grp)), c(grp$groups[["LandWebReport"]], "Mixed"))
  expect_error(reportingColors(given[-1L], grp), "no colour for: Pine")
})

test_that("the rebuilt maps keep the <rep>/vegTypeMap_year<YYYY>.tif path nrvtools parses", {
  withr::local_package("data.table")
  grp <- LandWebUtils::sppEquiv_groups(reportingEquiv(), "LandWeb", "LandWebReport")
  colors <- moduleFn(".reportingColors", new.env())(NULL, grp)
  buildMaps <- moduleFn(".buildReportingVegTypeMaps", new.env())

  repDir <- file.path(testPaths$outputPath, "reportingMaps", "mainSim", "rep01")
  dir.create(repDir, recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(reportingCohorts(), file.path(repDir, "cohortData_year0700.qs2"))
  terra::writeRaster(reportingPGM(), file.path(repDir, "pixelGroupMap_year0700.tif"),
                     datatype = "INT4U", overwrite = TRUE)
  vtm <- file.path(repDir, "vegTypeMap_year0700.tif") ## the species-level map: only its path is used

  outRoot <- file.path(testPaths$outputPath, "reportingMaps", "postprocess", "_reportingMaps")
  out <- buildMaps(vtm, outRoot, grp, colors, vegLeadingProportion = 0.75, mixedType = 2L)

  expect_identical(normalizePath(out), normalizePath(file.path(outRoot, "rep01", "vegTypeMap_year0700.tif")))
  expect_identical(cellLabels(terra::rast(out)), c("Bl_Spruce", "Decid", "Pine", "Mixed"))
  label <- paste0(basename(dirname(out)), ".", tools::file_path_sans_ext(basename(out)), "_ANC")
  expect_match(label, "^(.*?)\\.(.*)_year([0-9]+)_(.*)$") ## nrvtools' .parse_metric_labels() pattern

  ## without the saved cohorts there is nothing to rebuild from
  file.remove(file.path(repDir, "cohortData_year0700.qs2"))
  ## the message names both formats' paths
  expect_error(buildMaps(vtm, outRoot, grp, colors, 0.75, 2L),
               "missing:\n.*/rep01/cohortData_year0700.qs2 or .*/rep01/cohortData/year=0700/part-0.parquet$")
})

test_that("reporting-group maps are identical from .qs2 and parquet cohorts", {
  withr::local_package("data.table")
  grp <- LandWebUtils::sppEquiv_groups(reportingEquiv(), "LandWeb", "LandWebReport")
  colors <- moduleFn(".reportingColors", new.env())(NULL, grp)
  buildMaps <- moduleFn(".buildReportingVegTypeMaps", new.env())

  ## rep01 saved its cohorts as .qs2, rep02 as parquet (runs from 2.0.0.9046); the species as the simulation has them
  root <- file.path(testPaths$outputPath, "reportingMapsFormats")
  cd <- reportingCohorts()
  cd$speciesCode <- factor(cd$speciesCode)
  reps <- file.path(root, "mainSim", c("rep01", "rep02"))
  for (d in reps) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(reportingPGM(), file.path(d, "pixelGroupMap_year0700.tif"), datatype = "INT4U", overwrite = TRUE)
  }
  qs2::qs_save(cd, file.path(reps[1], "cohortData_year0700.qs2"))
  nrvtools::write_cohort_parquet(cd, file.path(reps[2], "cohortData"), 700)
  outRoot <- file.path(root, "postprocess", "_reportingMaps")
  out <- buildMaps(file.path(reps, "vegTypeMap_year0700.tif"), outRoot, grp, colors, 0.75, 2L)

  a <- terra::rast(out[1])
  b <- terra::rast(out[2])
  expect_identical(terra::values(b), terra::values(a))
  expect_identical(terra::levels(b), terra::levels(a))
  expect_identical(terra::coltab(b), terra::coltab(a))
  expect_identical(cellLabels(b), c("Bl_Spruce", "Decid", "Pine", "Mixed"))

  ## with both formats in a folder, the parquet cohorts are read: these .qs2 ones would make every pixel Decid
  decoy <- data.table::copy(cd)
  decoy$speciesCode <- factor(rep("Popu_tre", nrow(decoy)), levels = levels(cd$speciesCode))
  qs2::qs_save(decoy, file.path(reps[2], "cohortData_year0700.qs2"))
  out2 <- buildMaps(file.path(reps[2], "vegTypeMap_year0700.tif"), outRoot, grp, colors, 0.75, 2L)
  expect_identical(cellLabels(terra::rast(out2)), c("Bl_Spruce", "Decid", "Pine", "Mixed"))
  ## unless the cohort files InitMulti() found are given and the part file is not one of them (an outputsDF that
  ## lists the .qs2): the .qs2 is read
  out3 <- buildMaps(file.path(reps[2], "vegTypeMap_year0700.tif"), outRoot, grp, colors, 0.75, 2L,
                    cohortFiles = file.path(reps[2], "cohortData_year0700.qs2"))
  expect_identical(cellLabels(terra::rast(out3)), rep("Decid", 4L))
})

test_that("reporting-group maps read a .qs2 cohort file padded differently from the vegTypeMap's year", {
  withr::local_package("data.table")
  grp <- LandWebUtils::sppEquiv_groups(reportingEquiv(), "LandWeb", "LandWebReport")
  colors <- moduleFn(".reportingColors", new.env())(NULL, grp)
  buildMaps <- moduleFn(".buildReportingVegTypeMaps", new.env())
  root <- file.path(testPaths$outputPath, "reportingMapsPadded")
  cd <- reportingCohorts()
  cd$speciesCode <- factor(cd$speciesCode)
  d <- file.path(root, "mainSim", "rep01")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(reportingPGM(), file.path(d, "pixelGroupMap_year0700.tif"), datatype = "INT4U", overwrite = TRUE)
  padded <- file.path(d, "cohortData_year700.qs2")
  qs2::qs_save(cd, padded)
  out <- buildMaps(file.path(d, "vegTypeMap_year0700.tif"), file.path(root, "postprocess", "_reportingMaps"), grp,
                   colors, 0.75, 2L, cohortFiles = padded)
  expect_identical(cellLabels(terra::rast(out)), c("Bl_Spruce", "Decid", "Pine", "Mixed"))
})
