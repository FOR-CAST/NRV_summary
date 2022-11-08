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
  reqdPkgs = list("data.table", "dplyr", "fs", "future.apply", "ggplot2", "googledrive",
                  "landscapemetrics",
                  "PredictiveEcology/LandWebUtils@development",
                  "purrr", "raster", "sf",
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
    defineParameter("studyAreaNamesCol", "character", NA, NA, NA,
                    "column name used to identify names of subpolygons (features) the study area polygon."),
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
      #sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "HSI_PineMarten", "plot", .last())
      if (isTRUE(P(sim)$upload)) {
        sim <- scheduleEvent(sim, end(sim), "NRV_summary", "upload", .last())
      }
    },
    plot = {
      plotFun(sim) # example of a plotting function

      # ! ----- STOP EDITING ----- ! #
    },
    postprocess = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      sim <- landscapeMetrics(sim)
      #sim <- patchMetrics(sim) ## TODO: fix standAge maps @!!!!

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
    setdiff(P(sim)$timeSeriesTimes, mod$analysesOutputsTimes), padL = padL)), collapse = "|"),
                       mod$allouts, value = TRUE, invert = TRUE)

  ## TODO: inventory all files to ensure correct dir structure? compare against expected files?
  #filesUserHas <- fs::dir_ls(P(sim)$simOutputPath, recurse = TRUE, type = "file", glob = "*.qs")

  # filesNeeded <- data.table(file = mod$allouts2, exists = TRUE) ## TODO

  # if (!all(filesNeeded$exists)) {
  #   missing <- filesNeeded[exists == FALSE, ]$file
  #   stop("Some simulation files missing:\n", paste(missing, collapse = "\n"))
  # }

  stopifnot(length(mod$allouts2) == 2 * length(P(sim)$reps) * length(mod$analysesOutputsTimes))

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

## build landscape metrics tables from vegetation type maps (VTMs)
landscapeMetrics <- function(sim) {
  sA <- studyArea(sim$ml, 2)
  polyNames <- unique(sA[[P(sim)$studyAreaNamesCol]])

  .ncores <- pemisc::optimalClusterNum(5000, maxNumClusters = min(parallel::detectCores() / 2, 32L)) ## TODO: use module param
  options(future.availableCores.fallback = .ncores)

  vtmListByPoly <- future_lapply(mod$vtm, function(f) {
    byPoly <- lapply(polyNames, function(polyName) {
      poly <- sA[sA[[P(sim)$studyAreaNamesCol]] == polyName,]
      r <- raster::raster(f)
      rc <- raster::crop(r, poly)
      rcm <- raster::mask(r, poly)
      rcm
    })
    names(byPoly) <- paste(tools::file_path_sans_ext(basename(f)), polyNames , sep = "_") ## vegTypeMap_yearXXXX_polyName

    byPoly
  })
  names(vtmListByPoly) <- basename(dirname(mod$vtm)) ## repXX
  vtmListByPoly <- unlist(vtmListByPoly, recursive = FALSE, use.names = TRUE)

  labels <- purrr::transpose(strsplit(names(vtmListByPoly), "[.]"))
  labels1 <- unlist(labels[[1]])
  labels2 <- gsub("vegTypeMap_", "", unlist(labels[[2]]))
  labels2a <- purrr::transpose(strsplit(labels2, "_{1}"))
  labels2a1 <- unlist(labels2a[[1]])
  labels2a2 <- unlist(labels2a[[2]])

  vtmReps <- as.integer(gsub("rep", "", labels1))
  vtmTimes <- as.integer(gsub("year", "", labels2a1))
  vtmStudyAreas <- labels2a2

  funList <- list("lsm_l_area_mn",
                  "lsm_l_cohesion",
                  "lsm_l_condent",
                  "lsm_l_core_cv",
                  "lsm_l_ed",
                  "lsm_l_iji")
  names(funList) <- funList

  opt <- options(future.globals.maxSize = 5*1024^3) ## 5 GiB
  mod$fragStats <- future_lapply(funList, function(f) {
    fun <- get(f)

    frag_stat_df <- fun(vtmListByPoly) ## TODO: use non-default values?
    frag_stat_df <- mutate(frag_stat_df,
                           rep = vtmReps,
                           time = vtmTimes,
                           studyArea = vtmStudyAreas)
    frag_stat_df %>%
      group_by(time, studyArea) %>%
      summarise(N = length(value), mn = mean(value), sd = sd(value),
                se = sd / sqrt(N), ci = se * qt(0.975, N - 1))
  }, future.packages = "landscapemetrics")
  names(mod$fragStats) <- names(funList)
  options(opt)

  return(invisible(sim))
}

patchMetrics <- function(sim) {
  browser()
  ## TODO: identify problem with tsf maps -- all show tsf >> 600 years and same garbled values
  ## -- standAgeMaps all identical;
  ## -- ageMap completely garbled

  sA <- studyArea(sim$ml, 2)
  polyNames <- unique(sA[[P(sim)$studyAreaNamesCol]])

  .ncores <- pemisc::optimalClusterNum(5000, maxNumClusters = min(parallel::detectCores() / 2, 32L)) ## TODO: use module param
  options(future.availableCores.fallback = .ncores)

  tsfListByPoly <- future_lapply(mod$tsf, function(f) {
    byPoly <- lapply(polyNames, function(polyName) {
      poly <- sA[sA[[P(sim)$studyAreaNamesCol]] == polyName,]
      r <- raster::raster(f)
      rc <- raster::crop(r, poly)
      rcm <- raster::mask(r, poly)
      rcm
    })
    names(byPoly) <- paste(tools::file_path_sans_ext(basename(f)), polyNames , sep = "_") ## vegTypeMap_yearXXXX_polyName

    byPoly
  })
  names(tsfListByPoly) <- basename(dirname(mod$tsf)) ## repXX
  tsfListByPoly <- unlist(tsfListByPoly, recursive = FALSE, use.names = TRUE)

  labels <- purrr::transpose(strsplit(names(tsfListByPoly), "[.]"))
  labels1 <- unlist(labels[[1]])
  labels2 <- gsub("vegTypeMap_", "", unlist(labels[[2]]))
  labels2a <- purrr::transpose(strsplit(labels2, "_{1}"))
  labels2a1 <- unlist(labels2a[[1]])
  labels2a2 <- unlist(labels2a[[2]])

  tsfReps <- as.integer(gsub("rep", "", labels1))
  tsfTimes <- as.integer(gsub("year", "", labels2a1))
  tsfStudyAreas <- labels2a2

  if (FALSE) {
    # TODO: very few fires in studyAreaReporting!!
    sims <- file.path(unique(dirname(tsf)), "mySimOut_1000.qs")
    tmp_sim <- loadSimList(sims[1])
    raster::plot(tmp_sim$burnMap)
    raster::plot(studyArea(ml, 2), add = TRUE)
  }

  ## TODO: rework below to summarize for all landscape units
  ## TODO: use parallel
  lrgPtchs <- lapply(P(sim)$reps, function(rep) {
    rbindlist(lapply(mod$analysesOutputsTimes, function(year) {
      ids_tsf <- which(grepl(sprintf("rep%02d/rstTimeSinceFire_year%04d", rep, year), mod$tsf))
      ids_vtm <- which(grepl(sprintf("rep%02d/vegTypeMap_year%04d", rep, year), mod$vtm))
      ldt <- LargePatches(mod$tsf[ids_tsf], mod$vtm[ids_vtm], poly = studyArea(sim$ml, 2),
                          labelColumn = "shinyLabel", id = rep,
                          P(sim)$ageClassCutOffs, P(sim)$ageClasses, P(sim)$sppEquivCol, sim$sppEquiv)
      ldt[, time := year]
    }))
  })
  lrgPtchs_dt <- rbindlist(lrgPtchs)
  mod$vegTypes <- sort(unique(lrgPtchs_dt$vegCover))

  mod$patchAges <- lapply(mod$vegTypes, function(type) {
    lrgPtchs_dt[vegCover == type, ] %>%
      group_by(time, ageClass) %>%
      summarise(N = length(sizeInHa), mn = mean(sizeInHa), sd = sd(sizeInHa),
                se = sd / sqrt(N), ci = se * qt(0.975, N - 1),
                .groups = "keep")
  })
  names(mod$patchAges) <- mod$vegTypes

  return(invisible(sim))
}

### plotting
plot_over_time <- function(summary_df, ylabel) {
  ## TODO: each study area on separate facet or plot
  ggplot(summary_df, aes(x = time, y = mn, col = studyArea)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), width = 0.5) +
    ylab(ylabel)
}

plot_ptch_ages <- function(summary_df) {
  ggplot(summary_df, aes(x = time, y = mn, col = studyArea)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), width = 0.5) +
    ylab("Mean patch size (ha)") +
    facet_wrap(~ageClass, ncol = 1)
}

plotFun <- function(sim) {
  lapply(names(mod$fragStats), function(f) {
    ## TODO: use Plots
    #Plots(mod$fragStats[[f]], fn = plot_over_time, ylabel = substr(f, 7, nchar(f))) ## ??
    gg1 <- plot_over_time(mod$fragStats[[f]], substr(f, 7, nchar(f)))
    ggsave(file.path(outputPath(sim), "figures", paste0(f, ".png")), gg1)

    gg2 <- gg1 + facet_wrap(~studyArea)
    ggsave(file.path(outputPath(sim), "figures", paste0(f, "_facet.png")), gg2)
  })

  lapply(names(mod$patchAges), function(type) {
    ## TODO: use Plots
    #Plots(mod$patchAges[[type]], fn = plot_ptch_ages) ## ??
    gg <- plot_ptch_ages(mod$patchAges[[type]])
    ggsave(file.path(outputPath(sim), "figures", paste0(type, ".png")), gg)
  })

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
