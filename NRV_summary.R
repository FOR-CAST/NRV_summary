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
  reqdPkgs = list("data.table", "dplyr", "fs", "future.apply", "future.callr",
                  "ggforce", "ggplot2", "googledrive",
                  "landscapemetrics",
                  "PredictiveEcology/LandWebUtils@development (>= 0.1.5)",
                  "raster", "sf", "sp",
                  "PredictiveEcology/SpaDES.core@development (>= 1.1.1)"),
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

  mod$allouts <- fs::dir_ls(outputPath(sim), regexp = "vegType|TimeSince", recurse = 1, type = "file") |>
    grep("gri|png|txt|xml", x = _, value = TRUE, invert = TRUE)
  mod$allouts2 <- grep(paste(paste0("year", paddedFloatToChar(
    setdiff(c(0, P(sim)$timeSeriesTimes), mod$analysesOutputsTimes), padL = padL)), collapse = "|"),
    mod$allouts, value = TRUE, invert = TRUE)

  filesUserHas <- mod$allouts2

  dirsExpected <- file.path(outputPath(sim), sprintf("rep%02d", P(sim)$reps))
  filesExpected <- as.character(sapply(dirsExpected, function(d) {
    c(
      file.path(d, sprintf("rstTimeSinceFire_year%04d.tif", mod$analysesOutputsTimes)),
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

  mod$tsf <- gsub(".*vegTypeMap.*", NA, mod$allouts2) |>
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), x = _, value = TRUE)
  mod$vtm <- gsub(".*TimeSinceFire.*", NA, mod$allouts2) |>
    grep(paste(mod$analysesOutputsTimes, collapse = "|"), x = _, value = TRUE)

  mod$tsfTimeSeries <- gsub(".*vegTypeMap.*", NA, mod$allouts) |>
    grep(paste(P(sim)$timeSeriesTimes, collapse = "|"), x = _, value = TRUE)
  mod$vtmTimeSeries <- gsub(".*TimeSinceFire.*", NA, mod$allouts) |>
    grep(paste(P(sim)$timeSeriesTimes, collapse = "|"), x = _, value = TRUE)

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

calculateLandscapeMetrics <- function(summaryPolys, polyCol, vtm) {
  if (!is(summaryPolys, "sf"))
    summaryPolys <- sf::st_as_sf(summaryPolys)

  polyNames <- unique(summaryPolys[[polyCol]])

  funList <- list("lsm_l_area_mn",
                  "lsm_l_cohesion",
                  "lsm_l_condent",
                  "lsm_l_core_cv",
                  "lsm_l_ed",
                  "lsm_l_iji")
  names(funList) <- funList

  oldPlan <- plan(tweak(plan(), workers = pemisc::optimalClusterNum(5000, length(vtm))))
  on.exit(plan(oldPlan), add = TRUE)

  fragStats <- future.apply::future_lapply(vtm, function(f) {
    r <- terra::rast(f)
    byPoly <- lapply(polyNames, function(polyName) {
      subpoly <- summaryPolys[summaryPolys[[polyCol]] == polyName, ]
      rc <- terra::crop(r, subpoly)
      rcm <- terra::mask(rc, subpoly)
      rcm

      out <- lapply(funList, function(fun) {
        fn <- get(fun)

        fn(rcm)
      })
      names(out) <- funList
      out
    })
    names(byPoly) <- paste(tools::file_path_sans_ext(basename(f)), polyNames , sep = "_") ## vegTypeMap_yearXXXX_polyName

    byPoly
  }, future.packages = c("landscapemetrics", "sf", "terra"))
  names(fragStats) <- basename(dirname(vtm)) ## repXX

  fragStats <- purrr::transpose(lapply(fragStats, purrr::transpose)) ## puts fun names as outer list elements

  stopifnot(all(funList == names(fragStats)))

  frag_stat_df <- lapply(fragStats, function(x) {
    x <- unlist(x, recursive = FALSE, use.names = TRUE)

    labels <- purrr::transpose(strsplit(names(x), "[.]"))
    labels1 <- unlist(labels[[1]])
    labels2 <- gsub("vegTypeMap", "", unlist(labels[[2]]))
    labels2a <- purrr::transpose(strsplit(labels2, "_{1}"))
    labels2a2 <- unlist(labels2a[[2]]) ## year
    labels2a3 <- unlist(labels2a[[3]]) ## subpoly

    vtmReps <- as.integer(gsub("rep", "", labels1))
    vtmTimes <- as.integer(gsub("year", "", labels2a2))
    vtmStudyAreas <- labels2a3

    df <- do.call(rbind, x) |>
      mutate(rep = vtmReps, time = vtmTimes, poly = vtmStudyAreas) |>
      group_by(time, poly) |>
      summarise(N = length(value), mm = min(value), mn = mean(value, na.rm = TRUE), mx = max(value),
                sd = sd(value), se = sd / sqrt(N), ci = se * qt(0.975, N - 1))
  })
  names(frag_stat_df) <- funList

  return(frag_stat_df)
}

## build landscape metrics tables from vegetation type maps (VTMs)
landscapeMetrics <- function(sim) {
  ## current conditions
  vtmCC <- Cache(vegTypeMapGenerator,
                 x = sim$speciesLayers,
                 vegLeadingProportion = P(sim)$vegLeadingProportion,
                 mixedType = 2,
                 sppEquiv = sim$sppEquiv,
                 sppEquivCol = P(sim)$sppEquivCol,
                 colors = sim$sppColorVect,
                 doAssertion = FALSE)
  fname1 <- file.path(outputPath(sim), "vegTypeMap_year0000.tif")
  writeRaster(vtmCC, fname1, datatype = "INT1U", overwrite = TRUE)

  ## apply analysis to each of the reporting polygons
  md <- sim$ml@metadata
  cols <- which(grepl("analysisGroup", colnames(md)))
  rowIDs <- which(md[, ..cols] == currentModule(sim), arr.ind = TRUE)[, "row"]
  mod$rptPolyNames <- md[["layerName"]][rowIDs]
  lapply(mod$rptPolyNames, function(p) {
    message(crayon::magenta("Calculating landscape metrics for", p, "..."))

    rptPoly <- sim$ml[[p]]

    if (is(rptPoly, "Spatial")) {
      rptPoly <- st_as_sf(rptPoly)
    } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
      rptPoly <- st_collection_extract(rptPoly, "POLYGON")
    }
    rptPoly <- st_crop(rptPoly, studyArea(sim$ml, 3)) ## ensure cropped to studyArea

    rptPolyCol <- md[layerName == p, ][["columnNameForLabels"]]
    refCode <- paste0("lm_", md[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    vtm <- mod$vtm
    fileInfo <- file.info(vtm)[, c("size", "mtime")]
    mod[[refCodeCC]] <- suppressWarnings({
      Cache(calculateLandscapeMetrics, summaryPolys = rptPoly, polyCol = rptPolyCol, vtm = fname1,
            .cacheExtra = file.info(fname1)[, c("size", "mtime")])
    })
    mod[[refCode]] <- Cache(calculateLandscapeMetrics, summaryPolys = rptPoly, polyCol = rptPolyCol, vtm = vtm,
                            .cacheExtra = fileInfo)
  })

  return(invisible(sim))
}

## calculate areas for each patch (per species)
patchAreas <- function(vtm) {
  areas <- landscapemetrics::lsm_p_area(vtm)
  areas <- areas[areas$class != 0, ] ## class 0 has no forested vegetation (e.g., recently disturbed)
  spp <- raster::levels(vtm)[[1]]
  sppNames <- spp[match(areas$class, spp[["ID"]]), ][["values"]]

  areas <- mutate(areas, class = sppNames)

  return(areas)
}

## calculate median time since fire for each patch (per species)
patchAges <- function(vtm, tsf) {
  ptchs <- landscapemetrics::get_patches(vtm)[[1]] ## identify patches for each species (class)
  ptchs$class_0 <- NULL ## class 0 has no forested vegetation (e.g., recently disturbed)
  spp <- raster::levels(vtm)[[1]]
  spp$class <- paste0("class_", spp[["ID"]])
  names(ptchs) <- spp[match(names(ptchs), spp[["class"]]), ][["values"]]

  df <- rbindlist(lapply(names(ptchs), function(p) {
    ids <- which(!is.na(ptchs[[p]][]))
    data.frame(layer = 1L, level = "patch", class = p, id = ptchs[[p]][ids], metric = "tsf_mdn", tsf = tsf[ids]) |>
      group_by(layer, level, class, id, metric) |>
      summarise(value = median(tsf, na.rm = TRUE))
  }))

  return(df)
}

patchStats <- function(vtm, tsf, polyNames, summaryPolys, polyCol, funList) {
  t <- raster::raster(tsf)
  v <- raster::raster(vtm)
  byPoly <- lapply(polyNames, function(polyName) {
    message(paste("  vtm:", basename(vtm), "\n",
                  "  tsf:", basename(tsf)))
    subpoly <- summaryPolys[summaryPolys[[polyCol]] == polyName, ]

    tc <- raster::crop(t, subpoly)
    tcm <- raster::mask(tc, subpoly)
    tcm

    vc <- raster::crop(v, subpoly)
    vcm <- raster::mask(vc, subpoly)
    vcm

    out <- lapply(funList, function(fun) {
      message(paste("    ... running", fun, "for", polyName))

      fn <- get(fun)

      if (fun %in% c("patchAges")) {
        fn(vcm, tcm)
      } else {
        fn(vcm)
      }
    })
    names(out) <- funList
    out
  })
  names(byPoly) <- paste(tools::file_path_sans_ext(basename(vtm)), polyNames , sep = "_") ## vegTypeMap_yearXXXX_polyName

  byPoly
}

calculatePatchMetrics <- function(summaryPolys, polyCol, vtm, tsf) {
  if (!is(summaryPolys, "sf"))
    summaryPolys <- sf::st_as_sf(summaryPolys)

  polyNames <- unique(summaryPolys[[polyCol]])

  funList <- list("patchAges", "patchAreas")
  names(funList) <- funList

  oldPlan <- plan(tweak(plan(), workers = pemisc::optimalClusterNum(5000, length(vtm))))
  on.exit(plan(oldPlan), add = TRUE)

  ptch_stats <- future.apply::future_mapply(
    patchStats, vtm = vtm, tsf = tsf,
    MoreArgs = list(
     polyCol = polyCol,
     polyNames = polyNames,
     summaryPolys = summaryPolys,
     funList = funList
   ),
   SIMPLIFY = FALSE,
   future.globals = funList,
   future.packages = c("dplyr", "landscapemetrics", "raster", "sf") ## "terra"
  )
  names(ptch_stats) <- basename(dirname(vtm)) ## repXX

  ptch_stats <- purrr::transpose(lapply(ptch_stats, purrr::transpose)) ## puts fun names as outer list elements

  stopifnot(all(funList == names(ptch_stats)))

  ptch_stat_df <- lapply(ptch_stats, function(x) {
    x <- unlist(x, recursive = FALSE, use.names = TRUE)

    labels <- purrr::transpose(strsplit(names(x), "[.]"))
    labels1 <- unlist(labels[[1]])
    labels2 <- gsub("vegTypeMap", "", unlist(labels[[2]]))
    labels2a <- purrr::transpose(strsplit(labels2, "_{1}"))
    labels2a2 <- unlist(labels2a[[2]]) ## year
    labels2a3 <- unlist(labels2a[[3]]) ## subpoly

    vtmReps <- as.integer(gsub("rep", "", labels1))
    vtmTimes <- as.integer(gsub("year", "", labels2a2))
    vtmStudyAreas <- labels2a3

    df <- do.call(rbind, lapply(seq_along(x), function(i) {
      if (nrow(x[[i]]) == 0) {
        x[[i]] <- data.frame(layer = integer(0), level = character(0), class = character(0),
                             id = integer(0), metric = character(0), value = numeric(0))

      }
      mutate(x[[i]], rep = vtmReps[i], time = vtmTimes[i], poly = vtmStudyAreas[i]) |>
        group_by(class, time, poly, metric) |>
        summarise(
          N = length(value),
          mm = min(value, na.rm = TRUE),
          mn = mean(value, na.rm = TRUE),
          mx = max(value, na.rm = TRUE),
          sd = sd(value),
          se = sd / sqrt(N),
          ci = se * qt(0.975, N - 1)
        )
    }))
  })
  names(ptch_stat_df) <- funList

  return(ptch_stat_df)
}

patchMetrics <- function(sim) {
  ## current conditions
  vtmCC <- Cache(vegTypeMapGenerator,
                 x = sim$speciesLayers,
                 vegLeadingProportion = P(sim)$vegLeadingProportion,
                 mixedType = 2,
                 sppEquiv = sim$sppEquiv,
                 sppEquivCol = P(sim)$sppEquivCol,
                 colors = sim$sppColorVect,
                 doAssertion = FALSE)
  fname1 <- file.path(outputPath(sim), "vegTypeMap_year0000.tif")
  writeRaster(vtmCC, fname1, datatype = "INT1U", overwrite = TRUE)

  tsfCC <- sim$ml[["CC TSF"]]
  if (is(tsfCC, "PackedSpatRaster")) {
    tsfCC <- unwrap(tsfCC) ## TODO: why is this necessary???
  }
  fname2 <- file.path(outputPath(sim), "rstTimeSinceFire_year0000.tif")
  writeRaster(tsfCC, fname2, datatype = "INT1U", overwrite = TRUE)

  ## apply analysis to each of the reporting polygons
  md <- sim$ml@metadata
  cols <- which(grepl("analysisGroup", colnames(md)))
  rowIDs <- which(md[, ..cols] == currentModule(sim), arr.ind = TRUE)[, "row"]
  mod$rptPolyNames <- md[["layerName"]][rowIDs]
  lapply(mod$rptPolyNames, function(p) {
    message(crayon::magenta("Calculating patch metrics for", p, "..."))

    rptPoly <- sim$ml[[p]]

    if (is(rptPoly, "Spatial")) {
      rptPoly <- st_as_sf(rptPoly)
    } else if (is(rptPoly, "sf") && st_geometry_type(rptPoly, by_geometry = FALSE) != "POLYGON") {
      rptPoly <- st_collection_extract(rptPoly, "POLYGON")
    }
    rptPoly <- st_crop(rptPoly, studyArea(sim$ml, 3)) ## ensure cropped to studyArea
    rptPolyCol <- md[layerName == p, ][["columnNameForLabels"]]
    refCode <- paste0("pm_", md[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    ## CC
    fileInfo <- file.info(fname1, fname2)[, c("size", "mtime")]
    mod[[refCodeCC]] <- Cache(calculatePatchMetrics, tsf = fname2, vtm = fname1,
                              summaryPoly = rptPoly, polyCol = rptPolyCol,
                              .cacheExtra = fileInfo)

    ## simulation results
    fileInfo <- file.info(mod$tsf, mod$vtm)[, c("size", "mtime")]
    mod[[refCode]] <- Cache(calculatePatchMetrics, tsf = mod$tsf, vtm = mod$vtm,
                            summaryPoly = rptPoly, polyCol = rptPolyCol,
                            .cacheExtra = fileInfo)

    return(invisible(NULL))
  })

  return(invisible(sim))
}

### plotting
plot_over_time <- function(summary_df, ylabel, page = 1) {
  ggplot(summary_df, aes(x = time, y = mn)) +
    facet_wrap_paginate(~poly, ncol = 4, nrow = 3, page = page) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), width = 0.5) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab(ylabel)
}

plot_by_species <- function(summary_df, type = c("box", "violin"), page = 1) {
  ggplot(summary_df, aes(x = class, y = mn)) +
    facet_wrap_paginate(~poly, ncol = 4, nrow = 3, page = page) +
    switch(type,
           box = geom_boxplot(outlier.colour = "grey4", outlier.shape = 21, outlier.size = 1.0),
           violin = geom_violin(outlier.colour = "grey4")
    ) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 14)) +
    ylab(summary_df$metric)
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
      gg1 <- plot_over_time(mod[[refCode]][[f]], substr(f, 7, nchar(f))) +
        geom_hline(data = mod[[refCodeCC]][[f]], aes(yintercept = mn), col = "darkred", linetype = 2)
      nPages <- n_pages(gg1)
      lapply(seq_len(nPages), function(pg) {
        gg <- plot_over_time(mod[[refCode]][[f]], substr(f, 7, nchar(f)), page = pg) +
          geom_hline(data = mod[[refCodeCC]][[f]], aes(yintercept = mn), col = "darkred", linetype = 2)
        ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_p", pg, ".png")), gg,
               height = 10, width = 16)
      })
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
    refCode <- paste0("pm_", sim$ml@metadata[layerName == p, ][["shortName"]])
    refCodeCC <- paste0(refCode, "_CC")

    pngs2a <- lapply(names(mod[[refCode]]), function(f) {
      ## TODO: use Plots
      ggbox <- plot_by_species(mod[[refCode]][[f]], "box") +
        geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
      ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_box_plot.png")), ggbox,
             height = 10, width = 16)
    })

    pngs2b <- lapply(names(mod[[refCode]]), function(f) {
      ## TODO: use Plots
      ggvio <- plot_by_species(mod[[refCode]][[f]], "violin") +
        geom_point(data = mod[[refCodeCC]][[f]], col = "darkred", size = 2.5)
      ggsave(file.path(figurePath(sim), paste0(f, "_facet_by_", refCode, "_via_plot.png")), ggvio,
             height = 10, width = 16)
    })

    append(pngs2a, pngs2b)
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
  fsim <- file.path(outputPath(sim), paste0("simOutPreamble_", P(sim)$.studyAreaName, ".qs"))
  if (!file.exists(fsim)) {
    fsim <- file.path(tools::file_path_sans_ext(fsim), ".rds") ## fallback to rds if qs not used
  }
  tmp <- loadSimList(fsim)

  if (!suppliedElsewhere("ml", sim)) {
    sim$ml <- tmp$ml ## TODO: can't load ml objects from qs file !!
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    sim$sppEquiv <- tmp$sppEquiv
  }

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
