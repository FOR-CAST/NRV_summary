## A `postprocessEvents` value the module does not implement stops the run at init. Until 2.0.0.9027
## an unknown value was skipped without a word, and 'fd' was scheduled as a stub that stopped in
## browser() after the other summaries had run.

test_that("every implemented postprocess event is accepted, in any case", {
  check <- moduleFn(".checkPostprocessEvents", new.env())
  expect_no_error(check(c("am", "bc", "lm", "lw", "on", "pm")))
  expect_no_error(check(c("LM", "Pm")))
})

test_that("'fd' and unknown events stop, naming the values at fault", {
  check <- moduleFn(".checkPostprocessEvents", new.env())
  expect_error(check(c("lm", "fd")), "not implemented: 'fd'")
  expect_error(check(c("lm", "pn", "xx")), "not implemented: 'pn', 'xx'")
})

test_that("no module code calls browser()", {
  src <- file.path(moduleRoot, paste0(moduleName, ".R"))
  skip_if_not(file.exists(src), "module source not found")
  expect_false("browser" %in% all.names(parse(src, keep.source = FALSE)))
})
