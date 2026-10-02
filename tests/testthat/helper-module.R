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
