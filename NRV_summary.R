defineModule(sim, list(
  name = "NRV_summary",
  description = paste("NRV simulation post-processing and summary creation.",
                      "Produces 'X over time' summaries for multiple patch metrics."),
  keywords = c("NRV"),
  authors = c(
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut"))
  ),
  childModules = character(0),
  version = list(NRV_summary = "1.1.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "NRV_summary.Rmd"), ## same file
  reqdPkgs = list("data.table", "dplyr", "fs", "future.apply", "future.callr",
                  "ggforce", "ggplot2", "googledrive", "landscapemetrics",
                  "PredictiveEcology/LandR@development (>= 1.1.1)",
                  "PredictiveEcology/LandWebUtils@development (>= 0.1.5)",
                  "FOR-CAST/nrvtools (>= 0.0.21)",
                  "PredictiveEcology/pemisc@development (>= 0.0.4.9011)",
                  "raster", "sf", "sp",
                  "PredictiveEcology/SpaDES.core@development (>= 1.1.1)",
                  "terra"),
  parameters = bindrows(
    defineParameter("ageClasses", "character", LandWebUtils:::.ageClasses, NA, NA,
                    "descriptions/labels for age classes (seral stages)"),
    defineParameter("ageClassCutOffs", "integer", LandWebUtils:::.ageClassCutOffs, NA, NA,
                    "defines the age boundaries between age classes"),
    defineParameter("ageClassMaxAge", "integer", 400L, NA, NA,
                    "maximum possible age"),
    defineParameter("postprocessEvents", "character", c("lm", "pm"), NA, NA,
                    paste("Specify which subset of postprocessing events to run.",
                          "At least one of:",
                          "'lm' for default landscape metrics;",
                          "'pm' for default patch metrics;",
                          "'bc' for BC seral stage patch metrics;",
                          "'on' for ON patch metrics.")),
    defineParameter("reps", "integer", 1L:10L, 1L, NA_integer_,
                    paste("number of replicates/runs per study area.")),
    defineParameter("sieveThresh", "integer", 1L, NA_integer_, NA_integer_,
                    paste("threshold patch size (number of pixels) to use with `terra::sieve`",
                          "when creating seral stage maps")),
    defineParameter("sppEquivCol", "character", "EN_generic_short", NA, NA,
                    "The column in `sim$sppEquiv` data.table to use as a naming convention"),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "The column in `sim$sppEquiv` data.table to use as a naming convention"),
    defineParameter("summaryInterval", "integer", 100L, NA, NA,
                    "simulation time interval at which to take 'snapshots' used for summary analyses"),
    defineParameter("summaryPeriod", "integer", c(700L, 1000L), NA, NA,
                    "lower and upper end of the range of simulation times used for summary analyses"),
    defineParameter("timeSeriesTimes", "numeric", 601:650, NA, NA,
                    "simulation times for which to build time steries animations."),
    defineParameter("upload", "logical", FALSE, NA, NA,
                    "if TRUE, uses the `googledrive` package to upload figures."),
    defineParameter("uploadTo", "character", NA, NA, NA,
                    paste("if `upload = TRUE`, a Google Drive folder id corresponding to `.studyAreaName`.")),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0.0, 1.0,
                    "a number that defines whether a species is leading for a given pixel"),
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    paste("Human-readable name for the study area used - e.g., a hash of the study",
                          "area obtained using `reproducible::studyAreaName()`")),
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput("flammableMap", "SpatRaster",
                 desc = "binary flammability map. Required in single mode."),
    expectsInput("ml", "map",
                 desc = "map list object from preamble module (e.g., LandWeb_preamble)."),
    expectsInput("speciesLayers", "SpatRaster",
                 desc = "initial percent cover raster layers used for simulation."),
    expectsInput("sppColorVect", "character",
                 desc = paste("A named vector of colors to use for plotting.",
                              "The names must be in `sim$sppEquiv[[P(sim)$sppEquivCol]]`,",
                              "and should also contain a color for 'Mixed'")),
    expectsInput("sppEquiv", "data.table", NA, NA, NA,
                 desc = "table of species equivalencies. See `LandR::sppEquivalencies_CA`.")
  ),
  outputObjects = bindrows(
    createsOutput("ml", "map", "map list object"),
  )
))

## event types
#   - type `init` is required for initialization

doEvent.NRV_summary = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)

      if ("lm" %in% tolower(P(sim)$postprocessEvents)) {
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_lm", .last())
      }

      if ("pm" %in% tolower(P(sim)$postprocessEvents)) {
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_pm", .last())
      }

      if ("bc" %in% tolower(P(sim)$postprocessEvents)) {
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_bc", .last())
      }

      if ("on" %in% tolower(P(sim)$postprocessEvents)) {
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_on", .last())
      }

      sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess", .last())
      sim <- scheduleEvent(sim, end(sim), "NRV_summary", "plot", .last())

      if (isTRUE(P(sim)$upload)) {
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "upload", .last())
      }
    },
    plot = {
      plotFun(sim)
    },
    postprocess_lm = {
      sim <- landscapeMetrics(sim) ## TODO: warning: Number of classes must be >= 3, IJI = NA.
    },
    postprocess_pm = {
      sim <- patchMetrics(sim)
    },
    postprocess_bc = {
      sim <- makeSeralStageMapsBC(sim)
      sim <- patchMetricsSeralBC(sim)
    },
    postprocess_on = {
      ## TODO finalize implementation
      message("Ontario NRV metrics are not yet fully implemented.")
    },
    upload = {
      # ! ----- EDIT BELOW ----- ! #
      browser() ## TODO: split uploads based on P(sim)$postprocessEvent + test
      mod$files2upload <- set_names(mod$files2upload, basename(mod$files2upload))

      gid <- as_id(sim$uploadTo[[P(sim)$.studyAreaName]])
      prevUploaded <- drive_ls(gid)
      toUpload <- mod$files2upload[!(basename(mod$files2upload) %in% prevUploaded$name)]
      uploaded <- map(toUpload, ~ drive_upload(.x, path = gid))
      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

Init <- function(sim) {
  ## check for necessary output files -----------------------------------------------
  padL <- 4

  mod$analysesOutputsTimes <- analysesOutputsTimes(P(sim)$summaryPeriod, P(sim)$summaryInterval)

  cdpgm <- fs::dir_ls(outputPath(sim), regexp = "cohortData|pixelGroupMap", recurse = 1, type = "file") |>
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), x = _, value = TRUE)
  mod$allouts <- fs::dir_ls(outputPath(sim), regexp = "vegType|standAge", recurse = 1, type = "file") |>
    grep("gri|png|txt|xml", x = _, value = TRUE, invert = TRUE)
  mod$allouts2 <- grep(paste(paste0("year", paddedFloatToChar(
    setdiff(c(0, P(sim)$timeSeriesTimes), mod$analysesOutputsTimes), padL = padL)), collapse = "|"),
    mod$allouts, value = TRUE, invert = TRUE)

  filesUserHas <- c(cdpgm, mod$allouts2)

  dirsExpected <- file.path(outputPath(sim), sprintf("rep%02d", P(sim)$reps))
  filesExpected <- as.character(sapply(dirsExpected, function(d) {
    c(
      file.path(d, sprintf("cohortData_year%04d.qs", mod$analysesOutputsTimes)),
      file.path(d, sprintf("pixelGroupMap_year%04d.tif", mod$analysesOutputsTimes)),
      file.path(d, sprintf("standAgeMap_year%04d.tif", mod$analysesOutputsTimes)),
      file.path(d, sprintf("vegTypeMap_year%04d.tif", mod$analysesOutputsTimes))
    )
  }))

  filesNeeded <- data.frame(file = filesExpected, exists = filesExpected %in% filesUserHas)

  if (!all(filesNeeded$exists)) {
    missing <- filesNeeded[filesNeeded$exists == FALSE, ]$file
    stop(sum(!filesNeeded$exists), " simulation files appear to be missing:\n", paste(missing, collapse = "\n"))
  }

  mod$layerName <- gsub(mod$allouts2, pattern = paste0(".*", outputPath(sim)), replacement = "")
  mod$layerName <- gsub(mod$layerName, pattern = "[/\\]", replacement = "_")
  mod$layerName <- gsub(mod$layerName, pattern = "^_", replacement = "")

  mod$sam <- gsub(".*vegTypeMap.*", NA, mod$allouts2) |>
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), x = _, value = TRUE)
  mod$vtm <- gsub(".*standAgeMap.*", NA, mod$allouts2) |>
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), x = _, value = TRUE)

  mod$samTimeSeries <- gsub(".*vegTypeMap.*", NA, mod$allouts) |>
    grep(paste(P(sim)$timeSeriesTimes, collapse = "|"), x = _, value = TRUE)
  mod$vtmTimeSeries <- gsub(".*standAgeMap.*", NA, mod$allouts) |>
    grep(paste(P(sim)$timeSeriesTimes, collapse = "|"), x = _, value = TRUE)

  mod$flm <- file.path(outputPath(sim), "rstFlammable.tif")
  writeRaster(sim$flammableMap, mod$flm, overwrite = TRUE)

  ## cohortData and pixelGroupMap files
  mod$cd <- grep("cohortData", cdpgm, value = TRUE)
  mod$pgm <- grep("pixelGroupMap", cdpgm, value = TRUE)

  ## extract the reporting polygons to run the analyses on
  md <- sim$ml@metadata
  cols <- which(grepl("analysisGroup", colnames(md)))
  rowIDs <- which(md[, ..cols] == currentModule(sim), arr.ind = TRUE)[, "row"]
  mod$rptPolyNames <- md[["layerName"]][rowIDs]

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

## build landscape metrics tables from vegetation type maps (VTMs)
landscapeMetrics <- function(sim) {
  fvtm0 <- file.path(outputPath(sim), "vegTypeMap_year0000.tif")
  fvtm <- mod$vtm

  ## current conditions
  vtmCC <- Cache(vegTypeMapGenerator,
                 x = sim$speciesLayers,
                 vegLeadingProportion = P(sim)$vegLeadingProportion,
                 mixedType = 2,
                 sppEquiv = sim$sppEquiv,
                 sppEquivCol = P(sim)$sppEquivCol,
                 colors = sim$sppColorVect,
                 doAssertion = FALSE)
  writeRaster(vtmCC, fvtm0, datatype = "INT1U", overwrite = TRUE)

  md <- sim$ml@metadata

  ## sim$ml[[grep("\\(studyArea\\)", names(sim$ml), value = TRUE)]]
  studyArea3 <- map::studyArea(sim$ml, 3)

  funList <- default_landscape_metrics() ## TODO: pass this further up via parameter funList_lm

  oldPlan <- future::plan() |>
    tweak(workers = pemisc::optimalClusterNum(5000, length(fvtm))) |>
    future::plan()
  on.exit(future::plan(oldPlan), add = TRUE)

  lapply(mod$rptPolyNames, function(p, reportingPolygons, studyArea) {
    message(crayon::magenta("Calculating landscape metrics for", p, "..."))

    rptPoly <- reportingPolygons[[p]]

    if (is(rptPoly, "Spatial")) {
      rptPoly <- st_as_sf(rptPoly)
    } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
      rptPoly <- st_collection_extract(rptPoly, "POLYGON")
    }
    rptPoly <- st_crop(rptPoly, studyArea) ## ensure cropped to studyArea

    rptPolyCol <- md[layerName == p, ][["columnNameForLabels"]]
    refCode <- paste0("lm_", md[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    fileInfo <- file.info(fvtm0)[, c("size", "mtime")]
    mod[[refCodeCC]] <- suppressWarnings({
      Cache(calculateLandscapeMetrics, summaryPolys = rptPoly,
            polyCol = rptPolyCol, vtm = fvtm0, funList = funList,
            .cacheExtra = fileInfo)
    })
    lapply(names(mod[[refCodeCC]]), function(f) {
      write.csv(mod[[refCodeCC]][[f]], file.path(outputPath(sim), paste0(refCodeCC, "_", f, ".csv")), row.names = FALSE)
    })

    fileInfo <- file.info(vtm)[, c("size", "mtime")]
    mod[[refCode]] <- Cache(calculateLandscapeMetrics, summaryPolys = rptPoly,
                            polyCol = rptPolyCol, vtm = fvtm, funList = funList,
                            .cacheExtra = fileInfo)
    lapply(names(mod[[refCode]]), function(f) {
      write.csv(mod[[refCode]][[f]], file.path(outputPath(sim), paste0(refCode, "_", f, ".csv")), row.names = FALSE)
    })

    return(invisible(NULL))

  }, studyArea = studyArea3, reportingPolygons = sim$ml)

  return(invisible(sim))
}

patchMetrics <- function(sim) {
  fflm <- mod$flm
  fsam0 <- file.path(outputPath(sim), "standAgeMap_year0000.tif")
  fsam <- mod$sam
  fvtm0 <- file.path(outputPath(sim), "vegTypeMap_year0000.tif")
  fvtm <- mod$vtm

  ## current conditions
  vtmCC <- Cache(vegTypeMapGenerator,
                 x = sim$speciesLayers,
                 vegLeadingProportion = P(sim)$vegLeadingProportion,
                 mixedType = 2,
                 sppEquiv = sim$sppEquiv,
                 sppEquivCol = P(sim)$sppEquivCol,
                 colors = sim$sppColorVect,
                 doAssertion = FALSE)
  writeRaster(vtmCC, fvtm0, datatype = "INT1U", overwrite = TRUE)

  samCC <- if (is.null(sim$ml[["CC SAM"]])) sim$ml[["CC TSF"]] else sim$ml[["CC SAM"]]
  if (is(samCC, "PackedSpatRaster")) {
    samCC <- unwrap(samCC) ## TODO: why is this necessary? saveSimList wraps Spat* objects
  }
  writeRaster(samCC, fsam0, datatype = "INT1U", overwrite = TRUE)

  ## sim$ml[[grep("\\(studyArea\\)", names(sim$ml), value = TRUE)]]
  studyArea3 <- map::studyArea(sim$ml, 3)

  md <- sim$ml@metadata

  funList <- default_patch_metrics() ## TODO: pass this further up via parameter funList_pm

  oldPlan <- future::plan() |>
    tweak(workers = pemisc::optimalClusterNum(5000, length(fvtm))) |>
    future::plan()
  on.exit(future::plan(oldPlan), add = TRUE)

  lapply(mod$rptPolyNames, function(p, reportingPolygons, studyArea) {
    message(crayon::magenta("Calculating patch metrics for", p, "..."))

    rptPoly <- reportingPolygons[[p]]

    if (is(rptPoly, "Spatial")) {
      rptPoly <- st_as_sf(rptPoly)
    } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
      rptPoly <- st_collection_extract(rptPoly, "POLYGON")
    }
    rptPoly <- st_crop(rptPoly, studyArea) ## ensure cropped to studyArea
    rptPolyCol <- md[layerName == p, ][["columnNameForLabels"]]
    refCode <- paste0("pm_", md[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    ## CC
    fileInfo <- file.info(fvtm0, fsam0)[, c("size", "mtime")]
    dfl_cc <- Cache(calculatePatchMetrics, sam = fsam0, vtm = fvtm0, flm = fflm,
                    summaryPolys = rptPoly, polyCol = rptPolyCol,
                    funList = funList,
                    .cacheExtra = fileInfo)
    lapply(names(dfl_cc), function(f) {
      write.csv(dfl_cc[[f]], file.path(outputPath(sim), paste0(refCodeCC, "_", f, "_raw.csv")), row.names = FALSE)
    })

    mod[[refCodeCC]] <- summarizePatchMetrics(dfl_cc)
    lapply(names(mod[[refCodeCC]]), function(f) {
      write.csv(mod[[refCodeCC]][[f]], file.path(outputPath(sim), paste0(refCode, "_", f, ".csv")), row.names = FALSE)
    })

    ## simulation results
    fileInfo <- file.info(fsam, fvtm)[, c("size", "mtime")]
    dfl <- Cache(calculatePatchMetrics, sam = fsam, vtm = fvtm, flm = fflm,
                 summaryPolys = rptPoly, polyCol = rptPolyCol,
                 funList = funList,
                 .cacheExtra = fileInfo)
    lapply(names(dfl), function(f) {
      write.csv(dfl[[f]], file.path(outputPath(sim), paste0(refCode, "_", f, "_raw.csv")), row.names = FALSE)
    })

    mod[[refCode]] <- summarizePatchMetrics(dfl)
    lapply(names(mod[[refCode]]), function(f) {
      write.csv(mod[[refCode]][[f]], file.path(outputPath(sim), paste0(refCode, "_", f, ".csv")), row.names = FALSE)
    })

    return(invisible(NULL))
  }, studyArea = studyArea3, reportingPolygons = sim$ml)

  return(invisible(sim))
}

makeSeralStageMapsBC <- function(sim) {
  message(crayon::magenta("Creating seral stage maps ..."))

  ## sim$ml[[grep("\\(studyArea\\)", names(sim$ml), value = TRUE)]]
  studyArea3 <- map::studyArea(sim$ml, 3)
  NDTBEC <- sf::st_crop(sim$ml$`ecoregionLayer (NDTxBEC)`, studyArea3)
  fNDTBEC <- file.path(outputPath(sim), "NDTBEC.shp")
  sf::st_write(NDTBEC, fNDTBEC, append = FALSE, quiet = TRUE)
  rm(studyArea3, NDTBEC)

  fcd0 <- file.path(outputPath(sim), "rep01", "cohortData_year0000.qs")
  fpgm0 <- file.path(outputPath(sim), "rep01", "pixelGroupMap_year0000.tif")

  fcd <- c(fcd0, mod$cd)
  fpgm <- c(fpgm0, mod$pgm)

  oldPlan <- future::plan() |>
    tweak(workers = pemisc::optimalClusterNum(5000, length(fcd))) |>
    future::plan()
  on.exit(plan(oldPlan), add = TRUE)

  ssmFiles <- writeSeralStageMapBC(cd = fcd, pgm = fpgm, ndtbec = fNDTBEC)

  if (!is.na(P(sim)$sieveThresh)) {
    ## pass ssm through terra::sieve to merge singletons with neighbouring large patches?
    ## <https://rspatial.github.io/terra/reference/sieve.html>
    ssmFiles <- vapply(ssmFiles, function(f) {
      ## clumps < threshold merged with largest neighbour
      fs <- .suffix(f, sprintf("_sieve_%d", as.integer(P(sim)$sieveThresh))) ## TODO: use round() ??
      terra::sieve(rast(f), threshold = P(sim)$sieveThresh, filename = fs, overwrite = TRUE)
      fs
    }, character(1))
  }

  mod$ssm0 <- grep("seralStageMap_year0000.tif", ssmFiles, value = TRUE)
  mod$ssm <- grep("seralStageMap_year0000.tif", ssmFiles, invert = TRUE, value = TRUE)

  return(invisible(sim))
}

patchMetricsSeralBC <- function(sim) {
  fflm <- mod$flm
  fssm <- mod$ssm
  fssm0 <- mod$ssm0

  md <- sim$ml@metadata

  ## sim$ml[[grep("\\(studyArea\\)", names(sim$ml), value = TRUE)]]
  studyArea3 <- map::studyArea(sim$ml, 3)

  funList <- default_patch_metrics_seral() ## TODO: pass further up via parameter funList_bc

  rptPolygons <- sim$ml
  rptPolyNames <- mod$rptPolyNames

  oldPlan <- future::plan() |>
    tweak(workers = pemisc::optimalClusterNum(5000, length(fssm))) |>
    future::plan()
  on.exit(future::plan(oldPlan), add = TRUE)

  lapply(rptPolyNames, function(p, reportingPolygons, studyArea) {
    message(crayon::magenta("Calculating seral stage patch metrics for", p, "..."))

    rptPoly <- reportingPolygons[[p]]

    if (is(rptPoly, "Spatial")) {
      rptPoly <- st_as_sf(rptPoly)
    } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
      rptPoly <- st_collection_extract(rptPoly, "POLYGON")
    }
    rptPoly <- st_crop(rptPoly, studyArea) ## ensure cropped to studyArea
    rptPolyCol <- md[layerName == p, ][["columnNameForLabels"]]
    refCode <- paste0("sspm_", md[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    ## CC
    fileInfo <- file.info(fssm)[, c("size", "mtime")]
    dfl_cc <- Cache(calculatePatchMetricsSeral, ssm = fssm0, flm = fflm,
                    summaryPolys = rptPoly, polyCol = rptPolyCol,
                    funList = funList,
                    .cacheExtra = fileInfo)
    lapply(names(dfl_cc), function(f) {
      write.csv(dfl_cc[[f]], file.path(outputPath(sim), paste0(refCodeCC, "_", f, "_raw.csv")), row.names = FALSE)
    })
    mod[[refCodeCC]] <- summarizePatchMetricsSeral(dfl_cc)
    lapply(names(mod[[refCodeCC]]), function(f) {
      write.csv(mod[[refCodeCC]][[f]], file.path(outputPath(sim), paste0(refCodeCC, "_", f, ".csv")), row.names = FALSE)
    })

    ## simulation results
    fileInfo <- file.info(mod$ssm)[, c("size", "mtime")]
    dfl <- Cache(calculatePatchMetricsSeral, ssm = fssm, flm = fflm,
                 summaryPoly = rptPoly, polyCol = rptPolyCol,
                 funList = funList,
                 .cacheExtra = fileInfo)
    lapply(names(dfl), function(f) {
      write.csv(dfl[[f]], file.path(outputPath(sim), paste0(refCode, "_", f, "_raw.csv")), row.names = FALSE)
    })
    mod[[refCode]] <- summarizePatchMetricsSeral(dfl)
    lapply(names(mod[[refCode]]), function(f) {
      if (refCode == "sspm_NDTBEC" && f == "patchAreasSeral") {
        seral_table <- mod[[refCode]][[f]] |>
          na.omit() |>
          mutate(class = as.factor(class), poly = as.factor(poly),
                 mm = NULL, q1 = NULL, md = NULL, q3 = NULL, mx = NULL,
                 sd = NULL, cv = NULL, se = NULL, ci = NULL, n = NULL) |>
          ungroup() |>
          summarize(area = sum(N * mn, na.rm = TRUE), .by = c("class", "poly", "time")) |>
          mutate(totalArea = sum(area, na.rm = TRUE), .by = c("poly", "time")) |>
          summarize(
            minPctArea = 100 * min(area / totalArea, na.rm = TRUE),
            meanPctArea = 100 * mean(area / totalArea, na.rm = TRUE),
            maxPctArea = 100 * max(area / totalArea, na.rm = TRUE),
            .by = c("class", "poly")
          )

        write.csv(seral_table, file.path(outputPath(sim), "SeralTable.csv"), row.names = FALSE)
      }
      write.csv(mod[[refCode]][[f]], file.path(outputPath(sim), paste0(refCode, "_", f, ".csv")), row.names = FALSE)
    })

    return(invisible(NULL))
  }, studyArea = studyArea3, reportingPolygons = rptPolygons)

  return(invisible(sim))
}

### plotting
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  pngs_lm <- pngs_pm <- pngs_bc <- pngs_on <- list()

  if ("lm" %in% tolower(P(sim)$postprocessEvents)) {
    pngs_lm <- lapply(mod$rptPolyNames, function(p) {
      rptPoly <- sim$ml[[p]]

      if (is(rptPoly, "Spatial")) {
        rptPoly <- st_as_sf(rptPoly)
      } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
        rptPoly <- st_collection_extract(rptPoly, "POLYGON")
      }
      rptPolyCol <- sim$ml@metadata[layerName == p, ][["columnNameForLabels"]]
      refCode <- paste0("lm_", sim$ml@metadata[layerName == p, ][["shortName"]])
      refCodeCC <- paste0(refCode, "_CC")

      lapply(names(mod[[refCode]]), function(f) {
        ## TODO: use Plots
        gg1 <- plot_over_time(mod[[refCode]][[f]], substr(f, 7, nchar(f))) +
          geom_hline(data = mod[[refCodeCC]][[f]], aes(yintercept = mn), col = "darkred", linetype = 2)
        nPages <- n_pages(gg1)
        lapply(seq_len(nPages), function(pg) {
          gg <- plot_over_time(mod[[refCode]][[f]], substr(f, 7, nchar(f)), page = pg) +
            geom_hline(data = mod[[refCodeCC]][[f]], aes(yintercept = mn), col = "darkred", linetype = 2)
          ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_p", pg, ".png")), gg,
                 height = 10, width = 16)
        })
      }) |>
        unlist()
    })
  }

  if ("pm" %in% tolower(P(sim)$postprocessEvents)) {
    pngs_pm <- lapply(mod$rptPolyNames, function(p) {
      rptPoly <- sim$ml[[p]]

      if (is(rptPoly, "Spatial")) {
        rptPoly <- st_as_sf(rptPoly)
      } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
        rptPoly <- st_collection_extract(rptPoly, "POLYGON")
      }
      rptPolyCol <- sim$ml@metadata[layerName == p, ][["columnNameForLabels"]]
      refCode <- paste0("pm_", sim$ml@metadata[layerName == p, ][["shortName"]])
      refCodeCC <- paste0(refCode, "_CC")

      pngs_pm_a <- lapply(names(mod[[refCode]]), function(f) {
        ## TODO: use Plots
        ggbox1 <- plot_by_class(mod[[refCode]][[f]], "box") +
          geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
        nPages <- n_pages(ggbox1)
        lapply(seq_len(nPages), function(pg) {
          ggbox <- plot_by_class(mod[[refCode]][[f]], "box", page = pg) +
            geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
          ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_box_plot", "_p", pg, ".png")), ggbox,
                 height = 10, width = 16)
        })
      }) |>
        unlist()

      pngs_pm_b <- lapply(names(mod[[refCode]]), function(f) {
        ## TODO: use Plots
        ggvio1 <- plot_by_class(mod[[refCode]][[f]], "violin") +
          geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
        nPages <- n_pages(ggvio1)
        lapply(seq_len(nPages), function(pg) {
          ggvio <- plot_by_class(mod[[refCode]][[f]], "violin", page = pg) +
            geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
          ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_vio_plot", "_p", pg, ".png")), ggvio,
                 height = 10, width = 16)
        })
      }) |>
        unlist()

      c(pngs_pm_a, pngs_pm_b)
    })
  }

  if ("bc" %in% tolower(P(sim)$postprocessEvents)) {
    pngs_bc <- lapply(mod$rptPolyNames, function(p) {
      rptPoly <- sim$ml[[p]]

      if (is(rptPoly, "Spatial")) {
        rptPoly <- st_as_sf(rptPoly)
      } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
        rptPoly <- st_collection_extract(rptPoly, "POLYGON")
      }
      rptPolyCol <- sim$ml@metadata[layerName == p, ][["columnNameForLabels"]]
      refCode <- paste0("sspm_", sim$ml@metadata[layerName == p, ][["shortName"]])
      refCodeCC <- paste0(refCode, "_CC")

      pngs_bc_a <- lapply(names(mod[[refCode]]), function(f) {
        ## TODO: use Plots
        ggbox1 <- plot_by_class(mod[[refCode]][[f]], "box") +
          geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
        nPages <- n_pages(ggbox1)
        lapply(seq_len(nPages), function(pg) {
          ggbox <- plot_by_class(mod[[refCode]][[f]], "box", page = pg) +
            geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
          ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_box_plot", "_p", pg, ".png")), ggbox,
                 height = 10, width = 16)
        })
      }) |>
        unlist()

      pngs_bc_b <- lapply(names(mod[[refCode]]), function(f) {
        ## TODO: use Plots
        ggvio1 <- plot_by_class(mod[[refCode]][[f]], "violin") +
          geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
        nPages <- n_pages(ggvio1)
        lapply(seq_len(nPages), function(pg) {
          ggvio <- plot_by_class(mod[[refCode]][[f]], "violin", page = pg) +
            geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
          ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_vio_plot", "_p", pg, ".png")), ggvio,
                 height = 10, width = 16)
        })
      }) |>
        unlist()

      pngs_bc_c <- lapply(names(mod[[refCode]]), function(f) {
        ## TODO: use Plots
        gg1 <- plot_over_time_by_class(mod[[refCode]][[f]], f) +
          geom_hline(data = mod[[refCodeCC]][[f]], aes(yintercept = mn), linetype = 2)
        nPages <- n_pages(gg1)
        lapply(seq_len(nPages), function(pg) {
          gg <- plot_over_time_by_class(mod[[refCode]][[f]], f, page = pg) +
            geom_hline(data = mod[[refCodeCC]][[f]], aes(yintercept = mn), linetype = 2)
          ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_p", pg, ".png")), gg,
                 height = 10, width = 16)
        })
      }) |>
        unlist()

      c(pngs_bc_a, pngs_bc_b, pngs_bc_c)
    })
  }

  if ("on" %in% tolower(P(sim)$postprocessEvents)) {
    ## TODO
  }

  mod$files2upload <- c(
    unlist(pngs_lm),
    unlist(pngs_pm),
    unlist(pngs_bc),
    unlist(pngs_on)
  ) ## TODO: use registerOutputs()

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #
  fsim <- file.path(outputPath(sim), paste0("simOutPreamble_", P(sim)$.studyAreaName, ".qs"))
  if (!file.exists(fsim)) {
    fsim <- paste0(tools::file_path_sans_ext(fsim), ".rds") ## fallback to rds if qs not used
  }
  tmp <- loadSimList(fsim)

  if (!suppliedElsewhere("ml", sim)) {
    sim$ml <- tmp$ml ## TODO: can't load ml objects from qs file !!
  }

  if (!suppliedElsewhere("flammableMap", sim)) {
    sim$flammableMap <- tmp$flammableMap
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    sim$sppEquiv <- tmp$sppEquiv
  }

  ## TODO
  # if (!suppliedElsewhere("speciesLayers", sim)) {
  #   sim$speciesLayers <- tmp2$speciesLayers
  # }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

## older version of SpaDES.core used here doesn't have this function
if (packageVersion("SpaDES.core") < "2.0.2.9001") {
  figurePath <- function(sim) {
    file.path(outputPath(sim), "figures", current(sim)[["moduleName"]]) |>
      checkPath(create = TRUE)
  }
}
