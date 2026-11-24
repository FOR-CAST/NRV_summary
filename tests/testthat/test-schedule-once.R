## In mode = "single", `init` queues map_generators and save_single at the start, at the start of
## the summary period, at the end and at each timeSeriesTimes year, and each also reschedules itself
## summaryInterval later inside the summary period. Each distinct time must run once, and year 0
## keeps its place after every module's init and before time-0 growth (.first() + 1).

## The completed events of a mode = "single" run from 0 to 1000 with the shared fixture.
runSchedule <- function(timeSeriesTimes, out = file.path(testPaths$outputPath, "scheduleOnce")) {
  fx <- singleRunFixture()
  sim <- SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 1000),
    params = list(NRV_summary = list(
      mode = "single", summaryPeriod = c(700, 1000), summaryInterval = 50,
      timeSeriesTimes = timeSeriesTimes, ccStandAgeSource = "cohortData"
    )),
    modules = "NRV_summary",
    objects = fx$objects,
    paths = utils::modifyList(testPaths, list(outputPath = out))
  )
  cmp <- as.data.frame(SpaDES.core::completed(sim))
  cmp$eventTime <- as.numeric(cmp$eventTime)
  cmp
}

expectOncePerTime <- function(cmp, want) {
  for (type in c("map_generators", "save_single")) {
    times <- cmp$eventTime[cmp$eventType == type]
    expect_identical(sort(unique(times)), sort(as.numeric(want)), info = type)
    expect_equal(as.vector(table(times)), rep(1L, length(want)), info = type)
  }
}

test_that("each summary time has one map_generators and one save_single", {
  cmp <- runSchedule(timeSeriesTimes = 601:650)
  want <- c(0, 601:650, seq(700, 1000, 50))
  expectOncePerTime(cmp, want)

  y0 <- cmp[cmp$eventTime == 0 & cmp$eventType %in% c("map_generators", "save_single"), ]
  expect_equal(nrow(y0), 2L)
  expect_true(all(y0$eventPriority == SpaDES.core::.first() + 1))
})

test_that("a timeSeriesTimes year inside summaryPeriod does not cut the summary chain short", {
  ## 750 is both a timeSeriesTimes year and a summaryInterval step: init schedules it, so the
  ## reschedule from 700 skips it, and 750 must still reschedule 800 and the rest of the chain
  cmp <- runSchedule(timeSeriesTimes = c(601:650, 750))
  expectOncePerTime(cmp, c(0, 601:650, seq(700, 1000, 50)))
})

test_that("a timeSeriesTimes year off the summaryInterval grid is a one-off, not a second chain", {
  ## 725 is inside summaryPeriod but not on 700, 750, ...: it runs once and reschedules nothing,
  ## so 775, 825, ... never appear. The grid and the years before 700 are as without it.
  cmp <- runSchedule(timeSeriesTimes = c(601:650, 725))
  expectOncePerTime(cmp, c(0, 601:650, seq(700, 1000, 50), 725))
})

test_that("a saved year's vegTypeMap is INT1U, and saving raises no colour-table warning", {
  out <- withr::local_tempdir()
  warned <- character()
  withCallingHandlers(
    runSchedule(timeSeriesTimes = 601:602, out = out),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("color-table", warned)))
  f <- list.files(out, "^vegTypeMap_year[0-9]+\\.tif$", full.names = TRUE)
  expect_gt(length(f), 0L)
  for (x in f) {
    r <- terra::rast(x)
    expect_identical(terra::datatype(r), "INT1U", info = x)
    expect_true(terra::has.colors(r), info = x)
  }
})

test_that("a vegTypeMap code of 255 or more is refused: 255 is INT1U's NA value", {
  check <- moduleFn(".checkByteCodes", new.env())
  rtm <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1)
  expect_no_error(check(terra::rast(rtm, vals = c(1L, 5L, 254L))))
  expect_no_error(check(terra::rast(rtm, vals = c(1L, NA, 7L))))
  expect_error(check(terra::rast(rtm, vals = c(1L, 255L, 7L))), "below 255")
  expect_error(check(terra::rast(rtm, vals = c(1L, 300L, 7L))), "below 255")
})
