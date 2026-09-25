## The mode = "multi" outputs (parquet aggregates, envelope CSVs, animation, figures) are written
## outside the simList, so they reach a targets `_files` manifest only if registered with
## registerOutputs(). Until 2.0.0.9025 none were: the stage's `_files` target tracked 0 files.

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

test_that("postprocess registers the files the postprocess events wrote", {
  withr::local_package("SpaDES.core") ## as at run time: module code calls SpaDES.core unqualified
  mod <- new.env()
  f <- file.path(testPaths$outputPath, c("a.parquet", "b.csv"))
  file.create(f)
  mod$ppFiles <- c(f, f[1], file.path(testPaths$outputPath, "never-written.csv"))

  register <- moduleFn(".registerPostprocessOutputs", mod)
  sim <- register(SpaDES.core::simInit(paths = testPaths))

  registered <- normalizePath(SpaDES.core::outputs(sim)$file, winslash = "/", mustWork = FALSE)
  expect_setequal(registered, normalizePath(f, winslash = "/"))
  expect_length(mod$ppFiles, 0L) ## reset, so a second call does not register twice
})

test_that("the envelope CSVs record the files they write", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  mod$ppFiles <- character(0)
  writeCSVs <- moduleFn(".writeNrvSummaryCSVs", mod)
  sim <- SpaDES.core::simInit(paths = testPaths)
  env <- data.frame(metric = c("m1", "m2"), value = c(1, 2))

  writeCSVs(sim, env, "lm_ANSR", "ANSR")

  expect_length(mod$ppFiles, 3L) ## the combined CSV + one per metric
  expect_true(all(file.exists(mod$ppFiles)))
})

test_that("the plot event keeps the sim plotFun() returns", {
  ## plotFun() registers the PNGs on the sim it returns; calling it without assigning the result
  ## discarded every figure registration.
  doEvent <- get0(paste0("doEvent.", moduleName), mode = "function")
  body <- if (is.null(doEvent)) {
    readLines(file.path(moduleRoot, paste0(moduleName, ".R")))
  } else {
    deparse(doEvent)
  }
  expect_true(any(grepl("sim <- plotFun(sim)", body, fixed = TRUE)))
})
