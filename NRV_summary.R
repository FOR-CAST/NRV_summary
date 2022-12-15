defineModule(sim, list(
  name = "NRV_summary",
  description = paste("NRV simulation post-processing and summary creation.",
                      "Produces 'X over time' summaries for multiple patch metrics."),
  keywords = c("NRV"),
  authors = c(
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut"))
  ),
  childModules = character(0),
  version = list(NRV_summary = "0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "NRV_summary.Rmd"), ## same file
  reqdPkgs = list("data.table", "dplyr", "fs", "future.apply", "future.callr", "ggplot2", "googledrive",
                  "landscapemetrics",
                  "PredictiveEcology/LandWebUtils@development (>= 0.1.5)",
                  "raster", "sf", "sp",
                  "PredictiveEcology/SpaDES.core@development (>=1.1.0.9001)"),
  parameters = bindrows(
    defineParameter("ageClasses", "character", LandWebUtils:::.ageClasses, NA, NA, ## TODO: using 20 yr inc.
                    "descriptions/labels for age classes (seral stages)"),
    defineParameter("ageClassCutOffs", "integer", LandWebUtils:::.ageClassCutOffs, NA, NA,  ## TODO: using 20 yr inc.
                    "defines the age boundaries between age classes"),
    defineParameter("ageClassMaxAge", "integer", 400L, NA, NA,
                    "maximum possible age"),
    defineParameter("reps", "integer", 1L:10L, 1L, NA_integer_,
                    paste("number of replicates/runs per study area.")),
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
    expectsInput("ml", "map",
                 desc = "map list object from preamble module (e.g., LandWeb_preamble)."),
    expectsInput("speciesLayers", "RasterStack",
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
      sim <- scheduleEvent(sim, end(sim), "NRV_summary", "postprocess", .last())
      sim <- scheduleEvent(sim, end(sim), "NRV_summary", "plot", .last())

      if (isTRUE(P(sim)$upload)) {
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "upload", .last())
      }
    },
    plot = {
      plotFun(sim)
    },
    postprocess = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      sim <- landscapeMetrics(sim)
      sim <- patchMetrics(sim)

      # ! ----- STOP EDITING ----- ! #
    },
    upload = {
      # ! ----- EDIT BELOW ----- ! #
      browser() ## TODO
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

### template initialization
Init <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #

  padL <- 4

  mod$analysesOutputsTimes <- analysesOutputsTimes(P(sim)$summaryPeriod, P(sim)$summaryInterval)

  mod$allouts <- fs::dir_ls(outputPath(sim), regexp = "vegType|TimeSince", recurse = 1, type = "file") %>%
    grep("gri|png|txt|xml", ., value = TRUE, invert = TRUE)
  mod$allouts2 <- grep(paste(paste0("year", paddedFloatToChar(
    setdiff(c(0, P(sim)$timeSeriesTimes), mod$analysesOutputsTimes), padL = padL)), collapse = "|"),
    mod$allouts, value = TRUE, invert = TRUE)

  filesUserHas <- mod$allouts2

  dirsExpected <- file.path(outputPath(sim), sprintf("rep%02d", P(sim)$reps))
  filesExpected <- as.character(sapply(dirsExpected, function(d) {
    c(
      file.path(d, sprintf("rstTimeSinceFire_year%04d.tif", mod$analysesOutputsTimes)),
      file.path(d, sprintf("vegTypeMap_year%04d.grd", mod$analysesOutputsTimes))
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

  mod$tsf <- gsub(".*vegTypeMap.*", NA, mod$allouts2) %>%
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), ., value = TRUE)
  mod$vtm <- gsub(".*TimeSinceFire.*", NA, mod$allouts2) %>%
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), ., value = TRUE)

  mod$tsfTimeSeries <- gsub(".*vegTypeMap.*", NA, mod$allouts) %>%
    grep(paste(P(sim)$timeSeriesTimes, collapse = "|"), ., value = TRUE)
  mod$vtmTimeSeries <- gsub(".*TimeSinceFire.*", NA, mod$allouts) %>%
    grep(paste(P(sim)$timeSeriesTimes, collapse = "|"), ., value = TRUE)

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

calculateLandscapeMetrics <- function(summaryPolys, polyCol, vtm) {
  if (!is(summaryPolys, "sf"))
    summaryPolys <- sf::st_as_sf(summaryPolys)

  polyNames <- unique(summaryPolys[[polyCol]])

  ## vegetation type maps
  message("|_ loading vegetation type maps...")
  vtmListByPoly <- rasterListByPoly(files = vtm, poly = summaryPolys, names = polyNames,
                                    col = polyCol, filter = "vegTypeMap_") ## TODO:cache this
  vtmReps <- attr(vtmListByPoly, "reps")
  vtmTimes <- attr(vtmListByPoly, "times")
  vtmStudyAreas <- attr(vtmListByPoly, "polyNames")

  funList <- list("lsm_l_area_mn",
                  "lsm_l_cohesion",
                  "lsm_l_condent",
                  "lsm_l_core_cv",
                  "lsm_l_ed",
                  "lsm_l_iji")
  names(funList) <- funList

  fragStats <- future_lapply(funList, function(f) {
    fun <- get(f)

    frag_stat_df <- fun(vtmListByPoly) ## TODO: use non-default values?
    frag_stat_df <- mutate(frag_stat_df,
                           rep = vtmReps,
                           time = vtmTimes,
                           poly = vtmStudyAreas)
    frag_stat_df %>%
      group_by(time, poly) %>%
      summarise(N = length(value), mm = min(value), mn = mean(value), mx = max(value),
                sd = sd(value), se = sd / sqrt(N), ci = se * qt(0.975, N - 1))
  }, future.packages = "landscapemetrics")
  names(fragStats) <- names(funList)

  return(fragStats)
}

## build landscape metrics tables from vegetation type maps (VTMs)
landscapeMetrics <- function(sim) {
  .ncores <- pemisc::optimalClusterNum(5000, maxNumClusters = min(parallel::detectCores() / 2, 32L)) ## TODO: use module param
  opt <- options(future.availableCores.fallback = .ncores,
                 future.globals.maxSize = 0.26*length(P(sim)$reps)*1024^3) ## 50 reps needs ~13 GB
  on.exit({options(opt)}, add = TRUE)

  ## current conditions
  vtmCC <- Cache(vegTypeMapGenerator,
                 x = sim$speciesLayers,
                 vegLeadingProportion = P(sim)$vegLeadingProportion,
                 mixedType = 2,
                 sppEquiv = sim$sppEquiv,
                 sppEquivCol = P(sim)$sppEquivCol,
                 colors = sim$sppColorVect,
                 doAssertion = FALSE)
  fname1 <- file.path(outputPath(sim), "vegTypeMap_year0000.grd")
  raster::writeRaster(vtmCC, fname1, datatype = "INT1U", overwrite = TRUE)

  ## apply analysis to each of the reporting polygons
  md <- sim$ml@metadata
  rowIDs <- which(md == currentModule(sim), arr.ind = TRUE)[, "row"]
  mod$rptPolyNames <- md[["layerName"]][rowIDs]
  lapply(mod$rptPolyNames, function(p) {
    rptPoly <- sim$ml[[p]]

    if (is(rptPoly, "Spatial")) {
      rptPoly <- st_as_sf(rptPoly)
    } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
      rptPoly <- st_collection_extract(rptPoly, "POLYGON")
    }
    rptPolyCol <- md[layerName == p, ][["columnNameForLabels"]]
    refCode <- paste0("lm_", md[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    mod[[refCodeCC]] <- suppressWarnings({
      calculateLandscapeMetrics(summaryPolys = rptPoly, polyCol = rptPolyCol, vtm = fname1)
    })
    mod[[refCode]] <- calculateLandscapeMetrics(summaryPolys = rptPoly, polyCol = rptPolyCol, vtm = mod$vtm)
  })

  return(invisible(sim))
}

patchMetrics <- function(sim) {
  .ncores <- pemisc::optimalClusterNum(5000, maxNumClusters = min(parallel::detectCores() / 2, 32L)) ## TODO: use module param
  options(future.availableCores.fallback = .ncores)

  ## current conditions
  fname2 <- file.path(outputPath(sim), "rstTimeSinceFire_year0000.tif")
  fname1 <- paste0(tools::file_path_sans_ext(gsub("rstTimeSinceFire", "vegTypeMap", fname2)), ".grd")
  tsfCC <- sim$ml[["CC TSF"]]
  raster::writeRaster(tsfCC, fname2, datatype = "INT1U", overwrite = TRUE)

  ## apply analysis to each of the reporting polygons
  md <- sim$ml@metadata
  rowIDs <- which(md == currentModule(sim), arr.ind = TRUE)[, "row"]
  mod$rptPolyNames <- md[["layerName"]][rowIDs]
  foo <- lapply(mod$rptPolyNames, function(p) {
    rptPoly <- sim$ml[[p]]

    if (is(rptPoly, "Spatial")) {
      rptPoly <- st_as_sf(rptPoly)
    } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
      rptPoly <- st_collection_extract(rptPoly, "POLYGON")
    }
    rptPolyCol <- md[layerName == p, ][["columnNameForLabels"]]
    refCode <- paste0("patchAges_", md[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    ## -- BEGIN: to-rework
    ## TODO: rework below to summarize for all landscape units
    ## TODO: summarive current conditons
    ## TODO: use parallel

    subPolyNames <- unique(rptPoly[[rptPolyCol]])
    tsfListByPoly <- rasterListByPoly(files = mod$tsf, poly = rptPoly, names = subPolyNames,
                                      col = rptPolyCol, filter = "rstTimeSinceFire_")
    tsfReps <- attr(tsfListByPoly, "reps")
    tsfTimes <- attr(tsfListByPoly, "times")
    tsfStudyAreas <- attr(tsfListByPoly, "polyNames")
    tsfListByPoly <- lapply(tsfListByPoly, function(x) {
      x[] <- as.integer(pmin(P(sim)$ageClassMaxAge, x[]))
      x
    })

    ## CC
    sppEquivCol <- P(sim)$sppEquivCol
    age_veg <- CJ(ageClass = factor(P(sim)$ageClasses, levels = P(sim)$ageClasses), ## correct order
                  vegCover = unique(sim$sppEquiv[, ..sppEquivCol][[1]]))
    lrgPtchsCC_dt <- LargePatches(fname2, fname1, poly = sf::as_Spatial(rptPoly),
                                  labelColumn = rptPolyCol, id = rep,
                                  P(sim)$ageClassCutOffs, P(sim)$ageClasses, P(sim)$sppEquivCol, sim$sppEquiv)
    lrgPtchsCC_dt[, time := 0]
    lrgPtchsCC_dt[, ageClass := factor(ageClass, levels = P(sim)$ageClasses)] ## keep correct ageClass order
    vegTypesCC <- sort(unique(lrgPtchsCC_dt$vegCover))

    mod[[refCodeCC]] <- lapply(vegTypesCC, function(type) {
      lrgPtchsCC_dt[vegCover == type, ] %>%
        group_by(time, ageClass) %>%
        summarise(N = length(sizeInHa), mm = min(sizeInHa), mn = mean(sizeInHa), mx = max(sizeInHa),
                  sd = sd(sizeInHa), se = sd / sqrt(N), ci = se * qt(0.975, N - 1), .groups = "keep")
    })
    names(mod[[refCodeCC]]) <- vegTypesCC

    ## simulation results
    lrgPtchs <- lapply(P(sim)$reps, function(rep) {
      rbindlist(future_lapply(mod$analysesOutputsTimes, function(year) {
        ids_tsf <- which(grepl(sprintf("rep%02d/rstTimeSinceFire_year%04d", rep, year), mod$tsf))
        ids_vtm <- which(grepl(sprintf("rep%02d/vegTypeMap_year%04d", rep, year), mod$vtm))
        ldt <- LargePatches(mod$tsf[ids_tsf], mod$vtm[ids_vtm], poly = sf::as_Spatial(rptPoly),
                            labelColumn = rptPolyCol, id = rep,
                            P(sim)$ageClassCutOffs, P(sim)$ageClasses, P(sim)$sppEquivCol, sim$sppEquiv)
        ldt[, time := year]
      }, future.packages = c("LandWebUtils", "map"), future.seed = TRUE))
    }) ## TODO: currently quite slow...
    lrgPtchs_dt <- rbindlist(lrgPtchs)
    lrgPtchs_dt[, ageClass := factor(ageClass, levels = P(sim)$ageClasses)] ## keep correct ageClass order
    vegTypes <- sort(unique(lrgPtchs_dt$vegCover))

    mod[[refCode]] <- lapply(vegTypes, function(type) {
      lrgPtchs_dt[vegCover == type, ] %>%
        group_by(time, ageClass) %>%
        summarise(N = length(sizeInHa), mm = min(sizeInHa), mn = mean(sizeInHa), mx = max(sizeInHa),
                  sd = sd(sizeInHa), se = sd / sqrt(N), ci = se * qt(0.975, N - 1), .groups = "keep")
    })
    names(mod[[refCode]]) <- vegTypes

    ## -- END: to-rework
  })

  return(invisible(sim))
}

### plotting
plot_over_time <- function(summary_df, ylabel) {
  ggplot(summary_df, aes(x = time, y = mn)) +
    facet_wrap(~poly) +
    # geom_rect(aes(xmin = min(time), xmax = max(time), ymin = min(mm), ymax = max(mx)),
    #           fill = "grey", alpha = 0.1) + ## faceted plot already zoomed to range of data
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), width = 0.5) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab(ylabel)
}

plot_ptch_ages <- function(summary_df) {
  ggplot(summary_df, aes(x = time, y = mn)) +
    facet_wrap(~ageClass, ncol = 1) +
    # geom_rect(aes(xmin = min(time), xmax = max(time), ymin = min(mm), ymax = max(mx)),
    #           fill = "grey", alpha = 0.1) + ## faceted plot already zoomed to range of data
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), width = 0.5) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("Mean patch size (ha)")
}

plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  pngs1 <- lapply(mod$rptPolyNames, function(p) {
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
      gg <- plot_over_time(mod[[refCode]][[f]], substr(f, 7, nchar(f))) +
        geom_hline(data = mod[[refCodeCC]][[f]], aes(yintercept = mn), col = "red", linetype = 2)
      ggsave(file.path(outputPath(sim), "figures", paste0(f, "_facet_by_", refCode, ".png")), gg)
    })
  })

  pngs2 <- lapply(mod$rptPolyNames, function(p) {
    rptPoly <- sim$ml[[p]]

    if (is(rptPoly, "Spatial")) {
      rptPoly <- st_as_sf(rptPoly)
    } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
      rptPoly <- st_collection_extract(rptPoly, "POLYGON")
    }
    rptPolyCol <- sim$ml@metadata[layerName == p, ][["columnNameForLabels"]]
    refCode <- paste0("patchAges_", sim$ml@metadata[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    lapply(names(mod[[refCode]]),  function(spp) {
      ## TODO: use Plots
      gg <- plot_ptch_ages(mod[[refCode]][[spp]]) +
        geom_hline(data = mod[[refCodeCC]][[spp]], aes(yintercept = mn), col = "red", linetype = 2)
      ggsave(file.path(outputPath(sim), "figures", paste0(refCode, "_", spp, ".png")), gg)
    })
  })

  mod$files2upload <- c(unlist(pngs1), unlist(pngs2)) ## TODO: append these to sim outputs

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  tmp <- loadSimList(file.path(outputPath(sim), paste0("simOutPreamble_", P(sim)$.studyAreaName, ".qs")))

  if (!suppliedElsewhere("ml", sim)) {
    sim$ml <- tmp$ml ## TODO: can't load ml objects from qs file !!
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    sim$sppEquiv <- tmp$sppEquiv
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### add additional events as needed by copy/pasting from above
