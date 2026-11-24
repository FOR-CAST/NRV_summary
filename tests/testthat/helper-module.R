## Shared by the test files; testthat sources helper-*.R before them.

## A module function, with `mod` bound to `modEnv`. In CI the tests run against the package
## rendition (convertToPackage()), where the functions are in the namespace but `mod` is not; run
## by hand from the module, they are parsed out of <module>.R.
moduleFn <- function(name, modEnv) {
  fn <- get0(name, mode = "function")
  if (is.null(fn)) {
    src <- file.path(moduleRoot, paste0(moduleName, ".R"))
    skip_if_not(file.exists(src), "module source not found")
    defs <- new.env(parent = globalenv())
    for (x in parse(src, keep.source = FALSE)) {
      if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.call(x[[3L]]) &&
            identical(x[[3L]][[1L]], as.name("function"))) {
        eval(x, defs)
      }
    }
    fn <- get(name, envir = defs)
  }
  e <- new.env(parent = environment(fn))
  e$mod <- modEnv
  environment(fn) <- e
  fn
}

## The objects of a small mode = "single" run: two species on a 2 x 2 landscape, with the year-0
## cohorts. `objects` is ready for simInit(); `cohortData` is a fresh copy each call.
singleRunFixture <- function() {
  crs <- "EPSG:3400"
  rtm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = crs)
  spp <- c("Pice_gla", "Popu_tre")
  sppEquiv <- LandR::sppEquivalencies_CA
  sppEquiv <- sppEquiv[sppEquiv[["LandR"]] %in% spp, ]
  cohortData <- data.table::data.table(
    pixelGroup = c(1L, 1L, 2L, 3L),
    speciesCode = factor(c("Pice_gla", "Popu_tre", "Pice_gla", "Popu_tre"), levels = spp),
    age = c(80L, 40L, 150L, 20L),
    B = c(6000L, 2000L, 9000L, 1500L)
  )
  speciesLayers <- c(terra::rast(rtm, vals = c(60L, 60L, 100L, 0L)),
                     terra::rast(rtm, vals = c(40L, 40L, 0L, 100L)))
  names(speciesLayers) <- spp
  list(
    cohortData = cohortData,
    objects = list(
      cohortData = data.table::copy(cohortData),
      pixelGroupMap = terra::rast(rtm, vals = c(1L, 1L, 2L, 3L)),
      speciesLayers = speciesLayers,
      sppEquiv = sppEquiv,
      sppColorVect = LandR::sppColors(sppEquiv, "LandR", newVals = "Mixed", palette = "Accent"),
      studyAreaReporting = terra::as.polygons(terra::ext(rtm), crs = crs),
      flammableMap = terra::rast(rtm, vals = 1L)
    )
  )
}
