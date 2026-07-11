defineModule(sim, list(
  name = "NRV_summary",
  description = paste("NRV simulation post-processing and summary creation.",
                      "Produces summaries for multiple patch metrics and other indicators."),
  keywords = c("NRV"),
  authors = c(
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut"))
  ),
  childModules = character(0),
  version = list(NRV_summary = "2.0.0.9012"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  loadOrder = list(after = c("Biomass_core")),
  documentation = list("README.md", "NRV_summary.Rmd"), ## .md produced from .Rmd
  reqdPkgs = list(
    "data.table", "dplyr", "fs", "future.apply", "future.callr",
    "ggforce", "ggplot2", "gifski", "googledrive", "landscapemetrics", "qs2",
    "RColorBrewer", "sf", "terra", "tidyterra",
    "PredictiveEcology/LandR@development (>= 1.1.1)",
    "PredictiveEcology/LandWebUtils@development (>= 0.1.5)",
    "FOR-CAST/nrvtools (>= 0.2.0)",
    "PredictiveEcology/pemisc@development (>= 0.0.4.9011)",
    "PredictiveEcology/SpaDES.core@development (>= 3.0.3.9000)"
  ),
  parameters = bindrows(
    defineParameter("ageClasses", "character", LandWebUtils:::.ageClasses, NA, NA,
                    "descriptions/labels for age classes (seral stages)"),
    defineParameter("ageClassCutOffs", "integer", LandWebUtils:::.ageClassCutOffs, NA, NA,
                    "defines the age boundaries between age classes"),
    defineParameter("ageClassMaxAge", "integer", 400L, NA, NA,
                    "maximum possible age"),
    defineParameter("patchDirections", "integer", 4L, 4L, 8L,
                    paste("patch connectivity for the 'lw' (LandWeb summary) large-patch analysis,",
                          "passed to `nrvtools::largePatchCounts()` (via `landscapemetrics::get_patches`):",
                          "`4` = rook / 4-connected (matches v2's GDAL `polygonize`; the default),",
                          "`8` = queen / 8-connected. NOTE: v2 was fixed at 4-connectivity;",
                          "exposing 8 (queen) is a departure from v2.")),
    defineParameter("mixedType", "integer", 2L,
                    desc = paste("How to define mixed stands: `0L` for none; `1L` for any species admixture;",
                                 "`2L` for deciduous > conifer. See `LandR::vegTypeMapGenerator`.")),
    defineParameter("mode", "character", "single", NA, NA,
                    paste("use 'single' to run part of a simulation;",
                          "use 'multi' to run as part of postprocessing multiple runs.")),
    defineParameter("reuseAggregates", "logical", FALSE, NA, NA,
                    paste("(mode = 'multi') if TRUE, reuse a complete per-replicate `_aggregates`",
                          "parquet dataset instead of recomputing it, re-summarizing the surviving",
                          "parquets to regenerate the envelopes/CSVs/figures. Use to iterate on the",
                          "plots/CSVs without re-running the (~hours) landscape-metric aggregation;",
                          "leave FALSE for a fresh run.")),
    defineParameter("postprocessEvents", "character", c("lm", "pm"), NA, NA,
                    paste("Specify which subset of postprocessing events to run.",
                          "At least one of:",
                          "'am' for the stand-age time-series animation (GIF);",
                          "'bc' for BC seral stage patch metrics;",
                          "'fd' for forest degradation indicators;",
                          "'lm' for default landscape metrics;",
                          "'lw' for default LandWeb summaries",
                          "'pm' for default patch metrics;",
                          "'on' for ON patch metrics.")),
    defineParameter("reps", "integer", 1L:10L, 1L, NA_integer_,
                    paste("number of replicates/runs per study area.")),
    defineParameter("sieveThresh", "integer", 1L, NA_integer_, NA_integer_,
                    paste("threshold patch size (number of pixels) to use with `terra::sieve`",
                          "when creating seral stage maps")),
    defineParameter("simTimes", "numeric", c(NA, NA), NA, NA,
                    "Simulation start and end times when running in 'multi' mode."),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "The column in `sim$sppEquiv` data.table to use as a naming convention"),
    defineParameter("summaryInterval", "integer", 100L, NA, NA,
                    "simulation time interval at which to take 'snapshots' used for summary analyses."),
    defineParameter("summaryPeriod", "integer", start(sim) + c(700L, 1000L), NA, NA,
                    "lower and upper end of the range of simulation times used for summary analyses."),
    defineParameter("timeSeriesTimes", "numeric", start(sim) + 601:650, NA, NA,
                    "simulation times for which to build time steries animations."),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0.0, 1.0,
                    "a number that defines whether a species is leading for a given pixel"),
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by `SpaDES.core::Plots`, which can be optionally used here."),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    paste("Human-readable name for the study area used, or a hash of the study",
                          "area obtained using `reproducible::studyAreaName()`.")),
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput("cohortData", "data.table",
                 desc = "Required in single mode."),
    expectsInput("flammableMap", "SpatRaster",
                 desc = "binary flammability map (required with `type = 'single'`)"),
    expectsInput("pixelGroupMap", "SpatRaster",
                 desc = "Required in single mode."),
    expectsInput("reportingPolygons", "list",
                 desc = "reporting polygons for post-processing (required with `type = 'multi'`)"),
    expectsInput("speciesLayers", "SpatRaster",
                 desc = "initial percent cover raster layers used for simulation."),
    expectsInput("sppColorVect", "character",
                 desc = paste("A named vector of colors to use for plotting.",
                              "The names must be in `sim$sppEquiv[[P(sim)$sppEquivCol]]`,",
                              "and should also contain a color for 'Mixed'")),
    expectsInput("sppEquiv", "data.table", NA, NA, NA,
                 desc = "table of species equivalencies. See `LandR::sppEquivalencies_CA`."),
    expectsInput("studyAreaReporting", "SpatVector",
                 desc = "Required in single mode.")
  ),
  outputObjects = SpaDES.core:::._outputObjectsDF()
))

## event types
#   - type `init` is required for initialization

doEvent.NRV_summary = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      if (min(P(sim)$summaryPeriod) < start(sim) || max(P(sim)$summaryPeriod) > end(sim)) {
        stop("summaryPeriod values are outside the range of simulation times")
      }
      if (min(P(sim)$timeSeriesTimes) < start(sim) || max(P(sim)$timeSeriesTimes) > end(sim)) {
        stop("timeSeriesTimes values are outside the range of simulation times")
      }

      mod$analysesOutputsTimes <- analysesOutputsTimes(P(sim)$summaryPeriod, P(sim)$summaryInterval)

      if (P(sim)$mode == "single") {
        stopifnot(
          !is.null(sim$cohortData),
          !is.null(sim$pixelGroupMap),
          !is.null(sim$speciesLayers),
          !is.null(sim$sppColorVect),
          !is.null(sim$sppEquiv),
          !is.null(sim$studyAreaReporting)
        )

        sim <- scheduleEvent(sim, start(sim), "NRV_summary", "map_generators", .last())
        ## fmt: skip
        sim <- scheduleEvent(sim, P(sim)$summaryPeriod[1], "NRV_summary", "map_generators", .last())
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "map_generators", .last())

        sim <- scheduleEvent(sim, start(sim), "NRV_summary", "save_single", .last())
        sim <- scheduleEvent(sim, P(sim)$summaryPeriod[1], "NRV_summary", "save_single", .last())
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "save_single", .last())

        ## also generate + save the stand-age / veg-type maps at each timeSeriesTimes
        ## year, so the animation has its frames (read back in mode = "multi"). These
        ## typically fall outside summaryPeriod, so they are not covered by the
        ## summaryInterval reschedule; schedule map_generators before save_single so
        ## the maps exist when saved.
        for (tst in P(sim)$timeSeriesTimes) {
          sim <- scheduleEvent(sim, tst, "NRV_summary", "map_generators", .last())
          sim <- scheduleEvent(sim, tst, "NRV_summary", "save_single", .last())
        }
      } else if (P(sim)$mode == "multi") {
        stopifnot(!is.null(sim$reportingPolygons))

        ## reuse complete _aggregates parquet datasets (skip the ~2h landscape-metric recompute) so
        ## the postprocess events can iterate on plots/CSVs; read by .buildRepDataset() (process-local).
        options(NRV_summary.reuseAggregates = isTRUE(P(sim)$reuseAggregates))

        sim <- InitMulti(sim)

        if ("am" %in% tolower(P(sim)$postprocessEvents)) {
          sim <- scheduleEvent(sim, end(sim), "NRV_summary", "animation", .last())
        }

        if ("lm" %in% tolower(P(sim)$postprocessEvents)) {
          sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_lm", .last())
        }

        if ("pm" %in% tolower(P(sim)$postprocessEvents)) {
          sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_pm", .last())
        }

        if ("fd" %in% tolower(P(sim)$postprocessEvents)) {
          sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_fd", .last())
        }

        if ("lw" %in% tolower(P(sim)$postprocessEvents)) {
          sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_lw", .last())
        }

        if ("bc" %in% tolower(P(sim)$postprocessEvents)) {
          sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_bc", .last())
        }

        if ("on" %in% tolower(P(sim)$postprocessEvents)) {
          sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess_on", .last())
        }

        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess", .last())
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "plot", .last())
      }
    },
    map_generators = {
      mod$vegTypeMap <- LandR::vegTypeMapGenerator(
        sim$cohortData,
        sim$pixelGroupMap,
        P(sim)$vegLeadingProportion,
        mixedType = P(sim)$mixedType,
        sppEquiv = sim$sppEquiv,
        sppEquivCol = P(sim)$sppEquivCol,
        colors = sim$sppColorVect,
        doAssertion = getOption("LandR.assertions", TRUE)
      )

      mod$standAgeMap <- LandR::standAgeMapGenerator(
        sim$cohortData,
        sim$pixelGroupMap,
        weight = "biomass",
        doAssertion = getOption("LandR.assertions", TRUE)
      ) |>
        terra::mask(sim$studyAreaReporting)

      if (time(sim) >= P(sim)$summaryPeriod[1] && time(sim) < P(sim)$summaryPeriod[2]) {
        ## fmt: skip
        sim <- scheduleEvent(sim, time(sim) + P(sim)$summaryInterval, "NRV_summary", "map_generators", .last())
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
    postprocess_fd = {
      browser() ## TODO
    },
    postprocess_lw = {
      sim <- landWebMetrics(sim)
    },
    postprocess_bc = {
      sim <- makeSeralStageMapsBC(sim)
      sim <- patchMetricsSeralBC(sim)
    },
    postprocess_on = {
      ## TODO finalize implementation
      message("Ontario NRV metrics are not yet fully implemented.")
    },
    animation = {
      sim <- makeAnimation(sim)
    },
    save_single = {
      padYear <- paddedFloatToChar(time(sim), padL = ceiling(log10(end(sim) + 1)))

      ## objects to save at start of simulation ---------------------------------------------------
      if (time(sim) == start(sim)) {
        f_sppColorVect <- file.path(outputPath(sim), paste0("sppColorVect_year", padYear, ".qs2"))
        qs2::qs_save(sim$sppColorVect, f_sppColorVect)
        sim <- registerOutputs(f_sppColorVect, sim)

        f_sppEquiv <- file.path(outputPath(sim), paste0("sppEquiv_year", padYear, ".qs2"))
        qs2::qs_save(sim$sppEquiv, f_sppEquiv)
        sim <- registerOutputs(f_sppEquiv, sim)

        f_speciesLayers <- file.path(outputPath(sim), paste0("speciesLayers_year", padYear, ".tif"))
        terra::writeRaster(sim$speciesLayers, f_speciesLayers, datatype = "INT2U", overwrite = TRUE)
        sim <- registerOutputs(f_speciesLayers, sim)
      }

      ## objects to save during simulation --------------------------------------------------------
      ## stand-age + veg-type maps are saved at the summary times AND at each
      ## timeSeriesTimes year (the animation frames, read back in mode = "multi").
      ## cohortData + pixelGroupMap are only needed at the summary analysis times,
      ## so they are NOT written at the (many) timeSeriesTimes years.
      times_during <- c(start(sim), end(sim), mod$analysesOutputsTimes) |> unique() |> sort()
      times_maps <- c(times_during, P(sim)$timeSeriesTimes) |> unique() |> sort()

      if (time(sim) %in% times_maps) {
        f_standAgeMap <- file.path(outputPath(sim), paste0("standAgeMap_year", padYear, ".tif"))
        terra::writeRaster(mod$standAgeMap, f_standAgeMap, datatype = "INT2U", overwrite = TRUE)
        sim <- registerOutputs(f_standAgeMap, sim)

        f_vegTypeMap <- file.path(outputPath(sim), paste0("vegTypeMap_year", padYear, ".tif"))
        terra::writeRaster(mod$vegTypeMap, f_vegTypeMap, datatype = "INT2U", overwrite = TRUE)
        sim <- registerOutputs(f_vegTypeMap, sim)
      }

      if (time(sim) %in% times_during) {
        f_cohortData <- file.path(outputPath(sim), paste0("cohortData_year", padYear, ".qs2"))
        qs2::qs_save(sim$cohortData, f_cohortData)
        sim <- registerOutputs(f_cohortData, sim)

        f_pixelGroupMap <- file.path(outputPath(sim), paste0("pixelGroupMap_year", padYear, ".tif"))
        terra::writeRaster(sim$pixelGroupMap, f_pixelGroupMap, datatype = "INT4U", overwrite = TRUE)
        sim <- registerOutputs(f_pixelGroupMap, sim)

        if (time(sim) >= P(sim)$summaryPeriod[1] && time(sim) < P(sim)$summaryPeriod[2]) {
          ## fmt: skip
          sim <- scheduleEvent(sim, time(sim) + P(sim)$summaryInterval, "NRV_summary", "save_single", .last())
        }
      }

      ## objects to save at end of simulation -----------------------------------------------------
      if (time(sim) == end(sim)) {
        f_flammableMap <- file.path(outputPath(sim), paste0("flammableMap_year", padYear, ".tif"))
        terra::writeRaster(sim$flammableMap, f_flammableMap, datatype = "INT2U", overwrite = TRUE)
        sim <- registerOutputs(f_flammableMap, sim)
      }
    },
    noEventWarning(sim)
  )

  return(invisible(sim))
}

# event functions -------------------------------------------------------------------

InitMulti <- function(sim) {
  ## check for necessary output files -----------------------------------------------
  ## NOTE: don't load simLists -- slow and unreliable
  allReps <- sprintf("rep%02d", P(sim)$reps)
  padL <- ceiling(log10(P(sim)$simTimes[2] + 1))
  padYearStart <- paddedFloatToChar(P(sim)$simTimes[1], padL = padL)
  padYearEnd <- paddedFloatToChar(P(sim)$simTimes[2], padL = padL)

  ## all reps have same flammable map
  mod$flm <- file.path(outputPath(sim), allReps[1], paste0("flammableMap_year", padYearEnd, ".tif"))

  ## current-conditions reference = the sim's saved year-0 state (the deterministic
  ## initial condition, identical across reps -- read from rep 1). Read directly so
  ## the CC snapshot needs no regeneration from speciesLayers / no "CC SAM" input,
  ## and no write-before-read ordering between the landscape + patch metric events.
  mod$fvtm0 <- file.path(outputPath(sim), allReps[1], paste0("vegTypeMap_year", padYearStart, ".tif"))
  mod$fsam0 <- file.path(outputPath(sim), allReps[1], paste0("standAgeMap_year", padYearStart, ".tif"))
  ## current-conditions time-since-fire (burnSummaries output); age basis for the LandWeb summaries.
  mod$ftsf0 <- file.path(outputPath(sim), allReps[1], paste0("rstTimeSinceFire_year", padYearStart, ".tif"))

  cdpgm <- fs::dir_ls(
    outputPath(sim),
    regexp = "cohortData|pixelGroupMap",
    recurse = 1,
    type = "file"
  ) |>
    grep(paste0("(", paste0(allReps, collapse = "|"), ")"), x = _, value = TRUE) |>
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), x = _, value = TRUE)
  mod$allouts <- fs::dir_ls(
    outputPath(sim),
    regexp = "vegType|standAge",
    recurse = 1,
    type = "file"
  ) |>
    grep(paste0("(", paste0(allReps, collapse = "|"), ")"), x = _, value = TRUE) |>
    grep("gri|png|txt|xml", x = _, value = TRUE, invert = TRUE)
  mod$allouts2 <- paste(
    paste0(
      "year",
      paddedFloatToChar(
        setdiff(c(0, P(sim)$timeSeriesTimes), mod$analysesOutputsTimes),
        padL = padL
      )
    ),
    collapse = "|"
  ) |>
    grep(pattern = _, x = mod$allouts, value = TRUE, invert = TRUE)

  filesUserHas <- c(cdpgm, mod$allouts2)

  dirsExpected <- file.path(outputPath(sim), allReps)
  filesExpected <- as.character(sapply(dirsExpected, function(d) {
    c(
      file.path(d, sprintf("cohortData_year%04d.qs2", mod$analysesOutputsTimes)),
      file.path(d, sprintf("pixelGroupMap_year%04d.tif", mod$analysesOutputsTimes)),
      file.path(d, sprintf("standAgeMap_year%04d.tif", mod$analysesOutputsTimes)),
      file.path(d, sprintf("vegTypeMap_year%04d.tif", mod$analysesOutputsTimes))
    )
  }))

  filesNeeded <- data.frame(file = filesExpected, exists = filesExpected %in% filesUserHas)

  if (!all(filesNeeded$exists)) {
    missing <- filesNeeded[filesNeeded$exists == FALSE, ]$file
    stop(
      sum(!filesNeeded$exists),
      " simulation files appear to be missing:\n",
      paste(missing, collapse = "\n")
    )
  }

  mod$layerName <- gsub(mod$allouts2, pattern = paste0(".*", outputPath(sim)), replacement = "")
  mod$layerName <- gsub(mod$layerName, pattern = "[/\\]", replacement = "_")
  mod$layerName <- gsub(mod$layerName, pattern = "^_", replacement = "")

  mod$sam <- gsub(".*vegTypeMap.*", NA, mod$allouts2) |>
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), x = _, value = TRUE)
  mod$vtm <- gsub(".*standAgeMap.*", NA, mod$allouts2) |>
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), x = _, value = TRUE)

  ## time-since-fire per year (burnSummaries output), aligned with mod$vtm by rep/year path (see
  ## landWebMetrics(): the LandWeb summaries bin time-since-fire into age classes, matching v2).
  mod$tsf <- gsub("vegTypeMap", "rstTimeSinceFire", mod$vtm)

  mod$samTimeSeries <- gsub(".*vegTypeMap.*", NA, mod$allouts) |>
    grep(paste(P(sim)$timeSeriesTimes, collapse = "|"), x = _, value = TRUE)
  mod$vtmTimeSeries <- gsub(".*standAgeMap.*", NA, mod$allouts) |>
    grep(paste(P(sim)$timeSeriesTimes, collapse = "|"), x = _, value = TRUE)

  ## cohortData and pixelGroupMap files
  mod$cd <- grep("cohortData", cdpgm, value = TRUE)
  mod$pgm <- grep("pixelGroupMap", cdpgm, value = TRUE)

  ## extract the reporting polygons to run the analyses on
  mod$rptPolyNames <- names(sim$reportingPolygons)

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

## ---- arrow-native NRV summary helpers (nrvtools >= 0.1.0) --------------------------------------
## Replicated metrics are summarised through nrvtools' Arrow-native path (see the "Memory-bounded NRV
## summaries" vignette): each replicate's raw metric table is written to its own parquet partition
## under `<outputPath>/_aggregates/<refCode>/replicate=<rep>/`, then `summarize_nrv()` reduces across
## replicates by pushing the aggregation down to Arrow compute, so the per-replicate rows are never
## all held in memory at once. This replaces the former in-memory `calculateLandscapeMetrics()` /
## `summarizePatchMetrics()` reduction (both removed from nrvtools >= 0.1.0); the raw producers
## (`nrv_metrics_landscape()`, `calculatePatchMetrics()`, `calculatePatchMetricsSeral()`) now return
## the raw per-replicate long table (schema: level/class/metric/value/rep/time/poly) which
## `tidy_nrv_metrics()` binds and `summarize_nrv()` reduces to the mean/sd/min/max/... envelope.

## Post-processing outputs (parquet aggregates, figures, csv) live under
## `outputs/<studyArea>/postprocess/` -- a sibling of the mainSim rep dirs, which are still READ from
## `outputPath(sim)` (= outputs/<studyArea>/mainSim). `figures/` and `csv/` mirror the same
## `<kind>/<layer>/` sub-structure (kind = lm/pm/boxplots/histograms/...; layer = the full
## reporting-polygon-layer name), so a human can find a figure and its data side by side.
.ppRoot <- function(sim) {
  file.path(dirname(outputPath(sim)), "postprocess")
}

## `<postprocess>/_aggregates/<refCode>` -- the parquet dataset root for one refCode.
.nrvAggRoot <- function(sim, refCode) {
  file.path(.ppRoot(sim), "_aggregates", refCode)
}

## figure / csv output dir for one analysis `kind` and reporting `layer` (created on demand).
.ppFigDir <- function(sim, kind, layer, ...) {
  reproducible::checkPath(file.path(.ppRoot(sim), "figures", kind, layer, ...), create = TRUE)
}
.ppCsvDir <- function(sim, kind, layer, ...) {
  reproducible::checkPath(file.path(.ppRoot(sim), "csv", kind, layer, ...), create = TRUE)
}

## split a per-replicate map-file vector (`.../rep<NN>/<map>_year<YYYY>.tif`) into a named list by rep.
.filesByRep <- function(files) {
  split(files, basename(dirname(files)))
}

## TRUE iff `root` already holds a `replicate=<id>/*.parquet` partition for EVERY requested repID,
## i.e. the parquet dataset is complete and can be reused instead of recomputed.
.aggComplete <- function(root, repIDs) {
  dir.exists(root) &&
    all(vapply(repIDs, function(r) {
      length(list.files(file.path(root, paste0("replicate=", r)), pattern = "\\.parquet$")) > 0L
    }, logical(1L)))
}

## Build the parquet dataset for one refCode and return the across-replicate envelope.
## `compute_fn(repID)` returns the raw metric list for one replicate (from a raw nrvtools producer).
## `reuse = TRUE` skips the (expensive) recompute + rewrite when the dataset is already complete for
## `repIDs` (see `.aggComplete()`), re-summarizing the surviving parquets -- useful for iterating on
## the plots/CSVs without re-running the landscape-metric aggregation.
.buildRepDataset <- function(root, repIDs, compute_fn, studyArea = NULL, scenario = NULL,
                             id_cols = NULL, reuse = getOption("NRV_summary.reuseAggregates", FALSE)) {
  if (!(reuse && .aggComplete(root, repIDs))) {
    unlink(root, recursive = TRUE)
    for (repID in repIDs) {
      tidied <- tidy_nrv_metrics(compute_fn(repID), studyArea = studyArea, scenario = scenario)
      write_nrv_parquet(tidied, root, replicate = repID)
    }
  }
  ## id_cols = NULL -> summarize_nrv() default (per-time envelopes); the LandWeb summaries pass an
  ## explicit set excluding `time` so the NRV distribution pools across replicates AND summary years.
  summarize_nrv(root, id_cols = id_cols)
}

## Write the range-of-variation envelope for one refCode + reporting `layer` under csv/<kind>/<layer>/:
## a combined `<base>.csv` plus one `<base>_<metric>.csv` per metric. `kind` (lm/pm/sspm/lw) and the
## `_CC` suffix are recovered from `refCode`; the enclosing <kind>/<layer>/ dir supplies the context
## the flat filenames used to carry, and it mirrors the figure dir structure. The `lw` (LandWeb
## summary) envelope is split by analysis -- leadingProp -> csv/boxplots/<layer>/ (the leading
## boxplot data), the large-patch metrics -> csv/histograms/<layer>/.
.writeNrvSummaryCSVs <- function(sim, env, refCode, layer) {
  if (is.null(env) || !nrow(env)) {
    return(invisible(character(0)))
  }
  kind <- sub("_.*$", "", refCode)
  cc <- if (grepl("_CC$", refCode)) "_CC" else ""
  writeSet <- function(e, k, base) {
    if (is.null(e) || !nrow(e)) {
      return(invisible())
    }
    d <- .ppCsvDir(sim, k, layer)
    write.csv(e, file.path(d, paste0(base, cc, ".csv")), row.names = FALSE)
    for (m in unique(e$metric)) {
      write.csv(e[e$metric == m, ], file.path(d, paste0(base, cc, "_", m, ".csv")), row.names = FALSE)
    }
  }
  if (kind == "lw") {
    isLead <- env$metric == "leadingProp"
    writeSet(env[isLead, , drop = FALSE], "boxplots", "leading")
    writeSet(env[!isLead, , drop = FALSE], "histograms", "largePatches")
  } else {
    writeSet(env, kind, kind)
  }
  invisible()
}

## build landscape metric envelopes from vegetation type maps (VTMs)
landscapeMetrics <- function(sim) {
  fvtm0 <- mod$fvtm0 ## current-conditions VTM = the saved year-0 state (see InitMulti)
  fvtm <- mod$vtm
  studyAreaReporting <- sf::st_as_sf(sim$studyAreaReporting)

  funList <- default_landscape_metrics() ## TODO: pass this further up via parameter funList_lm

  oldPlan <- future::plan() |>
    tweak(workers = pemisc::optimalClusterNum(5000, length(fvtm))) |>
    future::plan()
  on.exit(future::plan(oldPlan), add = TRUE)

  vtmByRep <- .filesByRep(fvtm)

  lapply(
    mod$rptPolyNames,
    function(p, reportingPolygons, studyArea) {
      message(crayon::magenta("Calculating landscape metrics for", p, "..."))

      rptPoly <- reportingPolygons[[p]]

      if (is(rptPoly, "Spatial")) {
        rptPoly <- sf::st_as_sf(rptPoly)
      } else if (
        is(rptPoly, "sf") && sf::st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON"
      ) {
        rptPoly <- sf::st_collection_extract(rptPoly, "POLYGON")
      }
      rptPoly <- sf::st_crop(rptPoly, studyArea) ## ensure cropped to studyArea

      rptPolyCol <- "Name" ## label column set by LandWebUtils::buildReportingPolygons()
      refCode <- paste0("lm_", abbreviate(p, minlength = 8)) ## key output on the layer name (cf. bc event)
      refCodeCC <- paste0(refCode, "_CC")
      ## drop features with no grouping label: an NA polyName makes patch/landscape stats
      ## select an empty subpoly (`summaryPolys[[col]] == NA`) -> crop() NULL -> values(NULL).
      rptPoly <- rptPoly[!is.na(rptPoly[[rptPolyCol]]), ]
      if (nrow(rptPoly) == 0) {
        return(invisible(NULL)) ## no named features in this layer within the study area
      }

      ## raw per-replicate landscape metrics, Cached on the map file(s) it reads.
      lmRaw <- function(vtm) {
        Cache(
          nrv_metrics_landscape,
          summaryPolys = rptPoly,
          polyCol = rptPolyCol,
          vtm = vtm,
          funList = funList,
          .cacheExtra = file.info(vtm)[, c("size", "mtime")]
        )
      }

      ## current conditions: a single snapshot, summarised as one replicate.
      mod[[refCodeCC]] <- suppressWarnings(
        .buildRepDataset(
          .nrvAggRoot(sim, refCodeCC),
          repIDs = "CC",
          compute_fn = function(repID) lmRaw(fvtm0)
        )
      )

      ## simulation: one parquet partition per replicate, then reduce across reps.
      mod[[refCode]] <- .buildRepDataset(
        .nrvAggRoot(sim, refCode),
        repIDs = names(vtmByRep),
        compute_fn = function(repID) lmRaw(vtmByRep[[repID]])
      )

      .writeNrvSummaryCSVs(sim, mod[[refCode]], refCode, p)
      .writeNrvSummaryCSVs(sim, mod[[refCodeCC]], refCodeCC, p)

      return(invisible(NULL))
    },
    studyArea = studyAreaReporting,
    reportingPolygons = sim$reportingPolygons
  )

  return(invisible(sim))
}

patchMetrics <- function(sim) {
  fflm <- mod$flm
  ## current-conditions reference = the sim's saved year-0 state (see InitMulti):
  ## read directly, no regeneration from speciesLayers / no "CC SAM" input needed.
  fsam0 <- mod$fsam0
  fsam <- mod$sam
  fvtm0 <- mod$fvtm0
  fvtm <- mod$vtm

  studyAreaReporting <- sf::st_as_sf(sim$studyAreaReporting)
  funList <- default_patch_metrics() ## TODO: pass this further up via parameter funList_pm

  oldPlan <- future::plan() |>
    tweak(workers = pemisc::optimalClusterNum(5000, length(fvtm))) |>
    future::plan()
  on.exit(future::plan(oldPlan), add = TRUE)

  ## one parquet partition per replicate (the vtm/sam file vectors align by index).
  vtmByRep <- .filesByRep(fvtm)
  samByRep <- .filesByRep(fsam)

  lapply(
    mod$rptPolyNames,
    function(p, reportingPolygons, studyArea) {
      message(crayon::magenta("Calculating patch metrics for", p, "..."))

      rptPoly <- reportingPolygons[[p]]

      if (is(rptPoly, "Spatial")) {
        rptPoly <- st_as_sf(rptPoly)
      } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
        rptPoly <- st_collection_extract(rptPoly, "POLYGON")
      }
      rptPoly <- st_crop(rptPoly, studyArea) ## ensure cropped to studyArea
      rptPolyCol <- "Name" ## label column set by LandWebUtils::buildReportingPolygons()
      refCode <- paste0("pm_", abbreviate(p, minlength = 8)) ## key output on the layer name (cf. bc event)
      refCodeCC <- paste0(refCode, "_CC")
      ## drop features with no grouping label: an NA polyName makes patch/landscape stats
      ## select an empty subpoly (`summaryPolys[[col]] == NA`) -> crop() NULL -> values(NULL).
      rptPoly <- rptPoly[!is.na(rptPoly[[rptPolyCol]]), ]
      if (nrow(rptPoly) == 0) {
        return(invisible(NULL)) ## no named features in this layer within the study area
      }

      ## raw per-replicate patch metrics, Cached on the map files they read.
      pmRaw <- function(vtm, sam) {
        Cache(
          calculatePatchMetrics,
          sam = sam,
          vtm = vtm,
          flm = fflm,
          summaryPolys = rptPoly,
          polyCol = rptPolyCol,
          funList = funList,
          .cacheExtra = file.info(c(vtm, sam))[, c("size", "mtime")]
        )
      }

      ## current conditions
      mod[[refCodeCC]] <- .buildRepDataset(
        .nrvAggRoot(sim, refCodeCC),
        repIDs = "CC",
        compute_fn = function(repID) pmRaw(fvtm0, fsam0)
      )

      ## simulation results
      mod[[refCode]] <- .buildRepDataset(
        .nrvAggRoot(sim, refCode),
        repIDs = names(vtmByRep),
        compute_fn = function(repID) pmRaw(vtmByRep[[repID]], samByRep[[repID]])
      )

      .writeNrvSummaryCSVs(sim, mod[[refCode]], refCode, p)
      .writeNrvSummaryCSVs(sim, mod[[refCodeCC]], refCodeCC, p)

      return(invisible(NULL))
    },
    studyArea = studyAreaReporting,
    reportingPolygons = sim$reportingPolygons
  )

  return(invisible(sim))
}

## LandWeb summaries (ported v2 LandWeb_summary): leading-veg-by-age-class + large-patch counts.
## Age basis = time-since-fire (mod$tsf, from burnSummaries), matching v2. No flammable masking.
## The NRV distribution pools across replicates AND summary years, so summarize with `time` excluded
## from the id columns (idCols below); the plots (plotFun) are distributions, not time envelopes.
landWebMetrics <- function(sim) {
  ftsf0 <- mod$ftsf0
  ftsf <- mod$tsf
  fvtm0 <- mod$fvtm0
  fvtm <- mod$vtm

  studyAreaReporting <- sf::st_as_sf(sim$studyAreaReporting)
  funList <- default_landweb_metrics() ## TODO: pass this further up via parameter funList_lw
  idCols <- c("poly", "level", "class", "metric", "metric.1") ## pool across rep x summary year (no time)

  oldPlan <- future::plan() |>
    tweak(workers = pemisc::optimalClusterNum(5000, length(fvtm))) |>
    future::plan()
  on.exit(future::plan(oldPlan), add = TRUE)

  vtmByRep <- .filesByRep(fvtm)
  tsfByRep <- .filesByRep(ftsf)

  lapply(
    mod$rptPolyNames,
    function(p, reportingPolygons, studyArea) {
      message(crayon::magenta("Calculating LandWeb summaries for", p, "..."))

      rptPoly <- reportingPolygons[[p]]

      if (is(rptPoly, "Spatial")) {
        rptPoly <- sf::st_as_sf(rptPoly)
      } else if (
        is(rptPoly, "sf") && sf::st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON"
      ) {
        rptPoly <- sf::st_collection_extract(rptPoly, "POLYGON")
      }
      rptPoly <- sf::st_crop(rptPoly, studyArea) ## ensure cropped to studyArea
      rptPolyCol <- "Name" ## label column set by LandWebUtils::buildReportingPolygons()
      refCode <- paste0("lw_", abbreviate(p, minlength = 8)) ## key output on the layer name (cf. bc event)
      refCodeCC <- paste0(refCode, "_CC")
      rptPoly <- rptPoly[!is.na(rptPoly[[rptPolyCol]]), ]
      if (nrow(rptPoly) == 0) {
        return(invisible(NULL)) ## no named features in this layer within the study area
      }

      ## raw per-replicate LandWeb summaries, Cached on the map files they read.
      lwRaw <- function(vtm, tsf) {
        Cache(
          calculateLandWebMetrics,
          summaryPolys = rptPoly,
          polyCol = rptPolyCol,
          vtm = vtm,
          age = tsf,
          funList = funList,
          ageClassCutOffs = P(sim)$ageClassCutOffs,
          ageClasses = P(sim)$ageClasses,
          directions = P(sim)$patchDirections, ## 4 = rook (v2), 8 = queen; -> largePatchCounts()
          .cacheExtra = file.info(c(vtm, tsf))[, c("size", "mtime")]
        )
      }

      ## current conditions: a single (year-0) snapshot, treated as one replicate.
      mod[[refCodeCC]] <- .buildRepDataset(
        .nrvAggRoot(sim, refCodeCC),
        repIDs = "CC",
        compute_fn = function(repID) lwRaw(fvtm0, ftsf0),
        id_cols = idCols
      )

      ## simulation results: one parquet partition per replicate (vtm/tsf vectors align by index).
      mod[[refCode]] <- .buildRepDataset(
        .nrvAggRoot(sim, refCode),
        repIDs = names(vtmByRep),
        compute_fn = function(repID) lwRaw(vtmByRep[[repID]], tsfByRep[[repID]]),
        id_cols = idCols
      )

      .writeNrvSummaryCSVs(sim, mod[[refCode]], refCode, p)
      .writeNrvSummaryCSVs(sim, mod[[refCodeCC]], refCodeCC, p)

      return(invisible(NULL))
    },
    studyArea = studyAreaReporting,
    reportingPolygons = sim$reportingPolygons
  )

  return(invisible(sim))
}

makeSeralStageMapsBC <- function(sim) {
  message(crayon::magenta("Creating seral stage maps ..."))

  studyAreaReporting <- sf::st_as_sf(sim$studyAreaReporting)
  NDTBEC <- sim$reportingPolygons[["ecoregionLayer"]] |>
    sf::st_as_sf() |>
    sf::st_crop(studyAreaReporting)
  fNDTBEC <- file.path(outputPath(sim), "NDTBEC.shp")
  sf::st_write(NDTBEC, fNDTBEC, append = FALSE, quiet = TRUE)
  rm(studyAreaReporting, NDTBEC)

  fcd0 <- file.path(outputPath(sim), "rep01", "cohortData_year0000.qs2")
  fpgm0 <- file.path(outputPath(sim), "rep01", "pixelGroupMap_year0000.tif")

  fcd <- c(fcd0, mod$cd)
  fpgm <- c(fpgm0, mod$pgm)

  # oldPlan <- future::plan() |>
  #   tweak(workers = pemisc::optimalClusterNum(5000, length(fcd))) |>
  #   future::plan()
  oldPlan <- plan("callr", workers = pemisc::optimalClusterNum(5000, length(fcd)))
  on.exit(plan(oldPlan), add = TRUE)

  ssmFiles <- writeSeralStageMapBC(cd = fcd, pgm = fpgm, ndtbec = fNDTBEC)

  if (!is.na(P(sim)$sieveThresh)) {
    ## pass ssm through terra::sieve to merge singletons with neighbouring large patches?
    ## <https://rspatial.github.io/terra/reference/sieve.html>
    ssmFiles <- vapply(
      ssmFiles,
      function(f) {
        ## clumps < threshold merged with largest neighbour
        fs <- .suffix(f, sprintf("-sieve%d", as.integer(P(sim)$sieveThresh))) ## TODO: use round() ??
        terra::sieve(rast(f), threshold = P(sim)$sieveThresh, filename = fs, overwrite = TRUE)
        fs
      },
      character(1)
    )
  }

  mod$ssm0 <- grep("/seralStageMap_year0000.*[.]tif$", ssmFiles, value = TRUE)
  mod$ssm <- grep("/seralStageMap_year0000.*[.]tif$", ssmFiles, invert = TRUE, value = TRUE)

  return(invisible(sim))
}

patchMetricsSeralBC <- function(sim) {
  fflm <- mod$flm
  fssm0 <- mod$ssm0
  fssm <- mod$ssm
  studyAreaReporting <- sf::st_as_sf(sim$studyAreaReporting)

  funList <- default_patch_metrics_seral() ## TODO: pass further up via parameter funList_bc

  rptPolygons <- lapply(sim$reportingPolygons, sf::st_as_sf) ## converts to sf, keeping names
  rptPolyCols <- vapply(
    sim$reportingPolygons,
    FUN = attr,
    which = "field",
    FUN.VALUE = character(1)
  )
  rptPolyNames <- names(sim$reportingPolygons)

  oldPlan <- plan(workers = pemisc::optimalClusterNum(5000, length(fssm)))
  on.exit(future::plan(oldPlan), add = TRUE)

  ssmByRep <- .filesByRep(fssm)

  lapply(
    rptPolyNames,
    function(p, reportingPolygons, reportingPolygonCols, studyArea) {
      message(crayon::magenta("Calculating seral stage patch metrics for", p, "..."))

      rptPoly <- reportingPolygons[[p]]

      if (is(rptPoly, "sf") && sf::st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
        rptPoly <- sf::st_collection_extract(rptPoly, "POLYGON")
      }
      rptPoly <- sf::st_crop(rptPoly, studyArea) ## ensure cropped to studyArea
      rptPolyCol <- reportingPolygonCols[[p]]
      refCode <- paste0("sspm_", abbreviate(p, minlength = 8)) ## TODO: is this unique enough?
      refCodeCC <- paste0(refCode, "_CC")

      ## raw per-replicate seral patch metrics, Cached on the seral maps read.
      bcRaw <- function(ssm) {
        Cache(
          calculatePatchMetricsSeral,
          ssm = ssm,
          flm = fflm,
          summaryPolys = rptPoly,
          polyCol = rptPolyCol,
          funList = funList[[1]], ## TODO: temporarily, only patchAreasSeral
          .cacheExtra = file.info(ssm)[, c("size", "mtime")]
        )
      }

      ## current conditions
      mod[[refCodeCC]] <- .buildRepDataset(
        .nrvAggRoot(sim, refCodeCC),
        repIDs = "CC",
        compute_fn = function(repID) bcRaw(fssm0)
      )

      ## simulation results
      mod[[refCode]] <- .buildRepDataset(
        .nrvAggRoot(sim, refCode),
        repIDs = names(ssmByRep),
        compute_fn = function(repID) bcRaw(ssmByRep[[repID]])
      )

      .writeNrvSummaryCSVs(sim, mod[[refCode]], refCode, p)
      .writeNrvSummaryCSVs(sim, mod[[refCodeCC]], refCodeCC, p)

      ## SeralTable: range (min/mean/max over time) of each seral class's share of area, for the
      ## NDTxBEC reporting polygons. area = sum(n_reps * mean) reproduces the former sum(N * mn)
      ## (total patch area pooled across replicates); the replicate pooling cancels in the
      ## class/total proportion. `metric == "area"` is the landscapemetrics name for patchAreasSeral.
      if (refCode == "sspm_NDTBEC") {
        seral_table <- mod[[refCode]] |>
          dplyr::filter(.data$metric == "area") |>
          dplyr::select("class", "poly", "time", "n_reps", "mean") |>
          na.omit() |>
          dplyr::mutate(
            class = factor(.data$class, levels = seral_stages()),
            poly = as.factor(.data$poly)
          ) |>
          dplyr::summarize(
            area = sum(.data$n_reps * .data$mean, na.rm = TRUE),
            .by = c("class", "poly", "time")
          ) |>
          dplyr::mutate(totalArea = sum(.data$area, na.rm = TRUE), .by = c("poly", "time")) |>
          dplyr::summarize(
            minPctArea = 100 * min(.data$area / .data$totalArea, na.rm = TRUE),
            meanPctArea = 100 * mean(.data$area / .data$totalArea, na.rm = TRUE),
            maxPctArea = 100 * max(.data$area / .data$totalArea, na.rm = TRUE),
            .by = c("class", "poly")
          )

        write.csv(seral_table, file.path(outputPath(sim), "SeralTable.csv"), row.names = FALSE)
      }

      return(invisible(NULL))
    },
    studyArea = studyAreaReporting,
    reportingPolygons = rptPolygons,
    reportingPolygonCols = rptPolyCols
  )

  return(invisible(sim))
}

### stand-age time-series animation (ports the v2 LandWeb_summary `animation` event)
## Encoded with gifski (pure-Rust; no ImageMagick), so it avoids the ImageMagick
## cache-exhaustion that broke the v2 `animation::saveGIF` path -- see LandWeb#153.
## No system (policy.xml) configuration is required.
makeAnimation <- function(sim) {
  ## animate replicate 1's saved stand-age time series (the deterministic first rep,
  ## matching the CC snapshot); frames were masked to studyAreaReporting when saved.
  samFiles <- grep("rep01", mod$samTimeSeries, value = TRUE)
  if (length(samFiles) == 0L) {
    warning("NRV_summary animation: no rep01 standAgeMap time-series files found; skipping.")
    return(invisible(sim))
  }
  yrs <- as.integer(gsub(".*year0*([0-9]+)\\.tif$", "\\1", samFiles))
  ord <- order(yrs)
  samFiles <- samFiles[ord]
  yrs <- yrs[ord]

  ## age-class reclassification + colours (RdYlGn young -> old, matching v2's brewer.pal).
  ## `ageClassCutOffs` are the LOWER bound of each class (length == n); the final class
  ## runs to +Inf so old stands are never dropped.
  cutoffs <- P(sim)$ageClassCutOffs
  n <- length(cutoffs)
  ageClasses <- P(sim)$ageClasses[seq_len(n)]
  rcl <- cbind(cutoffs, c(cutoffs[-1], Inf), seq_len(n))
  pal <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(min(9L, max(3L, n)), "RdYlGn")
  )(n)
  names(pal) <- ageClasses

  ageClassFrame <- function(f) {
    r <- terra::classify(terra::rast(f), rcl, right = FALSE, include.lowest = TRUE)
    levels(r) <- data.frame(id = seq_len(n), ageClass = ageClasses)
    r
  }

  gifFile <- file.path(figurePath(sim), "standAge_animation.gif")
  gifski::save_gif(
    expr = {
      for (i in seq_along(samFiles)) {
        gg <- ggplot2::ggplot() +
          tidyterra::geom_spatraster(data = ageClassFrame(samFiles[i])) +
          ggplot2::scale_fill_manual(
            values = pal, na.value = "transparent", drop = FALSE, name = "age class"
          ) +
          ggplot2::labs(
            title = paste0(P(sim)$.studyAreaName, " — stand age"),
            subtitle = paste("year", yrs[i])
          ) +
          ggplot2::coord_sf(expand = FALSE) +
          ggplot2::theme_minimal()
        print(gg)
      }
    },
    gif_file = gifFile,
    width = 1200,
    height = 1200,
    delay = 1,
    progress = FALSE
  )

  message("NRV_summary: wrote stand-age animation (", length(samFiles), " frames) to ", gifFile)
  return(invisible(sim))
}

## LandWeb-summary figures for one reporting layer: the v2-form leading boxplots
## (figures/boxplots/<layer>/<subregion> <species>.png) and large-patch histograms
## (figures/histograms/<layer>/<size>/<subregion> <species>.png -- one file per species, four
## age-class panels), read from the raw per-replicate parquet with the current-condition overlay.
.saveLandWebFigs <- function(sim, p) {
  refCode <- paste0("lw_", abbreviate(p, minlength = 8))
  raw <- open_nrv_dataset(.nrvAggRoot(sim, refCode))
  if (is.null(raw)) {
    return(character(0))
  }
  raw <- as.data.frame(dplyr::collect(raw))
  if (!nrow(raw)) {
    return(character(0))
  }
  cc <- open_nrv_dataset(.nrvAggRoot(sim, paste0(refCode, "_CC")))
  cc <- if (!is.null(cc)) as.data.frame(dplyr::collect(cc)) else NULL

  ageClasses <- P(sim)$ageClasses
  saName <- P(sim)$.studyAreaName
  if (is.null(saName) || is.na(saName)) saName <- ""
  saPrefix <- if (nzchar(saName)) paste0(saName, " — ") else ""
  safe <- function(s) gsub("[/\\]", "-", s) ## filename-safe subregion / species
  ccFor <- function(poly, sp, met) {
    if (is.null(cc)) {
      NULL
    } else {
      cc[cc$poly == poly & cc$metric.1 == sp & cc$metric == met, , drop = FALSE]
    }
  }
  out <- character(0)

  ## Leading boxplots: one file per (subregion x species)
  lead <- raw[raw$metric == "leadingProp", , drop = FALSE]
  if (nrow(lead)) {
    dBox <- .ppFigDir(sim, "boxplots", p)
    for (poly in unique(lead$poly)) {
      for (sp in unique(lead$metric.1)) {
        d <- lead[lead$poly == poly & lead$metric.1 == sp, , drop = FALSE]
        if (!nrow(d) || all(d$value == 0, na.rm = TRUE)) next ## skip empty subregions / absent species
        gg <- nrvtools::plot_leading_boxplot(
          d,
          cc = ccFor(poly, sp, "leadingProp"),
          ageClasses = ageClasses,
          title = paste0(saPrefix, poly, " ", sp)
        )
        f <- file.path(dBox, paste0(safe(poly), " ", safe(sp), ".png"))
        ggsave(f, gg, width = 8, height = 6)
        out <- c(out, f)
      }
    }
  }

  ## Large-patch histograms: one file per (size x subregion x species), four age-class panels
  for (met in grep("^Npatch_ge", unique(raw$metric), value = TRUE)) {
    sz <- sub("^Npatch_ge(\\d+)ha$", "\\1", met)
    dSz <- .ppFigDir(sim, "histograms", p, sz)
    lp <- raw[raw$metric == met, , drop = FALSE]
    for (poly in unique(lp$poly)) {
      for (sp in unique(lp$metric.1)) {
        d <- lp[lp$poly == poly & lp$metric.1 == sp, , drop = FALSE]
        if (!nrow(d)) next
        gg <- nrvtools::plot_largepatch_histogram(
          d,
          cc = ccFor(poly, sp, met),
          ageClasses = ageClasses,
          xlab = paste("Number of patches greater than", sz, "ha"),
          title = paste0(saPrefix, poly, " ", sp, " (>=", sz, " ha)")
        )
        f <- file.path(dSz, paste0(safe(poly), " ", safe(sp), ".png"))
        ggsave(f, gg, width = 9, height = 7)
        out <- c(out, f)
      }
    }
  }
  out
}

### plotting
plotFun <- function(sim) {
  ## Envelope figures (lm/pm/bc) -> figures/<kind>/<layer>/{ribbon,boxplot}.png (faceted by the
  ## metric/class columns that vary). The LandWeb summaries (lw) are per-species boxplots / histograms
  ## via .saveLandWebFigs() -> figures/{boxplots,histograms}/<layer>/...
  saName <- P(sim)$.studyAreaName
  if (is.null(saName) || is.na(saName)) saName <- ""
  saPrefix <- if (nzchar(saName)) paste0(saName, " — ") else ""
  safe <- function(s) gsub("[/\\]", "-", s) ## filename-safe metric / subregion

  saveNrvPlots <- function(kind, p, ylab) {
    env <- mod[[paste0(kind, "_", abbreviate(p, minlength = 8))]]
    if (is.null(env) || !nrow(env)) {
      return(character(0))
    }
    d <- .ppFigDir(sim, kind, p)
    out <- character(0)
    ## One figure-set per metric: facet the (subregion x class) panels and paginate them
    ## across pages, so a large panel set becomes several PNGs (<metric>_<type>_p<pg>.png)
    ## instead of one crammed figure. Title carries the study area + metric name.
    for (met in unique(env$metric)) {
      sub <- env[env$metric == met, , drop = FALSE]
      if (!nrow(sub)) next
      ttl <- paste0(saPrefix, met)
      for (type in c("ribbon", "boxplot")) {
        gg1 <- plot_nrv_envelope(
          sub, type = type, facet = c("poly", "class", "metric.1"),
          ylab = ylab, title = ttl, page = 1
        )
        if (is.null(gg1)) next
        nPages <- tryCatch(ggforce::n_pages(gg1), error = function(e) 1L)
        if (is.null(nPages) || is.na(nPages)) nPages <- 1L
        for (pg in seq_len(nPages)) {
          gg <- plot_nrv_envelope(
            sub, type = type, facet = c("poly", "class", "metric.1"),
            ylab = ylab, title = ttl, page = pg
          )
          f <- file.path(d, paste0(safe(met), "_", type, "_p", pg, ".png"))
          ggsave(f, gg, height = 10, width = 16)
          out <- c(out, f)
        }
      }
    }
    out
  }

  events <- tolower(P(sim)$postprocessEvents)
  pngs <- character(0)

  if ("lm" %in% events) {
    pngs <- c(pngs, unlist(lapply(mod$rptPolyNames, function(p) {
      saveNrvPlots("lm", p, ylab = "landscape metric value")
    })))
  }
  if ("pm" %in% events) {
    pngs <- c(pngs, unlist(lapply(mod$rptPolyNames, function(p) {
      saveNrvPlots("pm", p, ylab = "patch metric value")
    })))
  }
  if ("lw" %in% events) {
    pngs <- c(pngs, unlist(lapply(mod$rptPolyNames, function(p) .saveLandWebFigs(sim, p))))
  }
  if ("bc" %in% events) {
    pngs <- c(pngs, unlist(lapply(mod$rptPolyNames, function(p) {
      saveNrvPlots("sspm", p, ylab = "seral patch area (ha)")
    })))
  }

  if (length(pngs)) {
    sim <- registerOutputs(pngs, sim)
  }
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  ## nothing here

  return(invisible(sim))
}
