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
  reqdPkgs = list("data.table", "dplyr", "fs", "ggplot2", "googledrive", "landscapemetrics",
                  "PredictiveEcology/LandWebUtils@development",
                  "raster", "sf",
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
    defineParameter("summaryInterval", "integer", 100L, NA, NA,
                    "simulation time interval at which to take 'snapshots' used for summary analyses"),
    defineParameter("summaryPeriod", "integer", c(700L, 1000L), NA, NA,
                    "lower and upper end of the range of simulation times used for summary analyses"),
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
  mod$allouts2 <- grep(paste(paste0("year", paddedFloatToChar(P(sim)$timeSeriesTimes, padL = padL)), collapse = "|"),
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
  browser()
  vtmReps <- as.integer(substr(basename(dirname(mod$vtm)), 4, 5)) ## keep as integer for calculations
  vtmTimes <- as.numeric(substr(basename(mod$vtm), 16, 19))
  vtmList <- lapply(mod$vtm, function(f) {
    r <- raster::raster(f)
    rc <- raster::crop(r, studyArea(sim$ml, 2))
    rcm <- raster::mask(r, studyArea(sim$ml, 2))
    rcm
  })
  funList <- list("lsm_l_area_mn",
                  "lsm_l_cohesion",
                  "lsm_l_condent",
                  "lsm_l_core_cv",
                  "lsm_l_ed",
                  "lsm_l_iji")
  names(funList) <- funList

  mod$fragStats <- lapply(funList, function(f) {
    fun <- get(f)

    frag_stat_df <- fun(vtmList) ## TODO: use non-default values?
    frag_stat_df <- mutate(frag_stat_df, rep = vtmReps, time = vtmTimes)
    frag_stat_summary_df <- frag_stat_df %>%
      group_by(time) %>%
      summarise(N = length(value), mn = mean(value), sd = sd(value),
                se = sd / sqrt(N), ci = se * qt(0.975, N - 1))
  })
  names(sim$fragStats) <- funList

  return(invisible(sim))
}

patchMetrics <- function(sim) {
  browser()
  ## TODO: identify problem with tsf maps -- all show tsf >> 600 years and same garbled values
  ## -- standAgeMaps all identical;
  ## -- ageMap completely garbled
  tsfReps <- as.integer(substr(basename(dirname(mod$tsf)), 4, 5)) ## keep as integer for calculations
  tsfTimes <- as.numeric(substr(basename(mod$tsf), 22, 25))
  tsfList <- lapply(tsf, function(f) {
    r <- raster::raster(f)
    rc <- raster::crop(r, studyArea(ml, 2))
    rcm <- raster::mask(r, studyArea(ml, 2))
    rcm
  })

  if (FALSE) {
    # TODO: very few fires in studyAreaReporting!!
    sims <- file.path(unique(dirname(tsf)), "mySimOut_1000.qs")
    tmp_sim <- loadSimList(sims[1])
    raster::plot(tmp_sim$burnMap)
    raster::plot(studyArea(ml, 2), add = TRUE)
  }

  allReps <- seq(startRep, length.out = nReps)
  lrgPtchs <- lapply(allReps, function(rep) {
    rbindlist(lapply(analysesOutputsTimes, function(year) {
      ids_tsf <- which(grepl(sprintf("rep%02d/rstTimeSinceFire_year%04d", rep, year), tsf))
      ids_vtm <- which(grepl(sprintf("rep%02d/vegTypeMap_year%04d", rep, year), vtm))
      ldt <- LargePatches(tsf[ids_tsf], vtm[ids_vtm], poly = studyArea(ml, 2), labelColumn = "shinyLabel",
                          id = rep, ageClassCutOffs, ageClasses, sppEquivCol, sppEquiv)
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
  ggplot(summary_df, aes(x = time, y = mn)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), width = 0.5) +
    ylab(ylabel)
}

plot_ptch_ages <- function(summary_df) {
  ggplot(summary_df, aes(x = time, y = mn)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), width = 0.5) +
    ylab("Mean patch size (ha)") +
    facet_wrap(~ageClass, ncol = 1)
}

plotFun <- function(sim) {
  lapply(names(sim$fragStats), function(f) {
    ## TODO: use Plots
    #Plots(sim$fragStats[[f]], fn = plot_over_time, ylabel = substr(f, 7, nchar(f))) ## ??
    gg <- plot_over_time(sim$fragStats[[f]], substr(f, 7, nchar(f)))
    ggsave(file.path(outputPath(sim), "figures", paste0(f, ".png")), gg)
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