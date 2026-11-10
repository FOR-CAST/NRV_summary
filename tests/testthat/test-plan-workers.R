## .planWithWorkers() sets the worker count on the caller's future backend for the metric events.
## Until 2.0.0.9028 the events re-tweaked the plan inline: under `sequential` that warned, and under a
## plan set as `plan(<backend>, workers = n)` it stopped, because future appends tweaks and so passed
## `workers` twice. The BC seral event called `plan(workers = n)`, which changes nothing.

test_that("a sequential plan is left alone, without a warning", {
  old <- future::plan(future::sequential)
  withr::defer(future::plan(old))
  planWithWorkers <- moduleFn(".planWithWorkers", new.env())

  expect_no_warning(planWithWorkers(3L))
  expect_s3_class(future::plan(), "sequential")
  expect_equal(future::nbrOfWorkers(), 1)
})

test_that("the worker count is set on a backend the caller left at its default", {
  skip_if_not(parallelly::supportsMulticore(), "multicore futures not supported here")
  old <- future::plan(future::multicore)
  withr::defer(future::plan(old))
  planWithWorkers <- moduleFn(".planWithWorkers", new.env())

  planWithWorkers(3L)
  expect_equal(future::nbrOfWorkers(), 3)
})

test_that("the worker count replaces one the caller set, and the previous plan comes back", {
  skip_if_not(parallelly::supportsMulticore(), "multicore futures not supported here")
  old <- future::plan(future::multicore, workers = 2L)
  withr::defer(future::plan(old))
  planWithWorkers <- moduleFn(".planWithWorkers", new.env())

  prev <- planWithWorkers(3L)
  expect_equal(future::nbrOfWorkers(), 3)
  future::plan(prev)
  expect_equal(future::nbrOfWorkers(), 2)
})
