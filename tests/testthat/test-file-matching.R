## InitMulti() finds each replicate's saved files by the year and replicate in the path's own fields.
## Unanchored greps matched a year such as "0" or "2820", or a replicate such as "rep1", anywhere in
## the path, so a folder name could select (or drop) files of other years and replicates.

test_that(".fileYear reads the year from the file name, NA where there is none", {
  fileYear <- moduleFn(".fileYear", new.env())
  f <- c("o/rep01/vegTypeMap_year0700.tif", "o/2820/standAgeMap_year0000.tif", "o/cohortData.qs2")
  expect_equal(fileYear(f), c(700L, 0L, NA))
  expect_type(fileYear(f), "integer")
})

test_that(".inYears ignores digits in folder names", {
  inYears <- moduleFn(".inYears", new.env())
  f <- "o_2820/rep01/vegTypeMap_year0700.tif"
  expect_false(inYears(f, 2820))
  expect_true(inYears(f, 700))
  expect_true(inYears(f, c(0, 700)))
  expect_false(inYears("o/cohortData.qs2", 0))
})

test_that(".inReps matches whole path components only", {
  inReps <- moduleFn(".inReps", new.env())
  f <- c("o/rep01/a.tif", "o/rep10/a.tif", "o/rep1/a.tif")
  expect_equal(inReps(f, "rep01"), c(TRUE, FALSE, FALSE))
  expect_false(inReps("o/rep10/a.tif", "rep1"))
  expect_equal(inReps(f, c("rep01", "rep10")), c(TRUE, TRUE, FALSE))
})

test_that("aux.xml sidecars carry no year", {
  fileYear <- moduleFn(".fileYear", new.env())
  inYears <- moduleFn(".inYears", new.env())
  f <- "o/rep01/pixelGroupMap_year0700.tif.aux.xml"
  expect_true(is.na(fileYear(f)))
  expect_false(inYears(f, 700))
})
