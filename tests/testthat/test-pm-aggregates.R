## The pm aggregates behind the patch metric envelopes and per-class pages (2.0.0.9032): absent
## classes get a class area of 0, and a patch-level metric's envelope spans the replicates' means.

## class areas (and the mean patch area) of the classes in `classes`, subregion A, one replicate
new_ca_replicate <- function(rep, classes, time) {
  d <- expand.grid(class = classes, time = time, stringsAsFactors = FALSE)
  rbind(
    data.frame(level = "class", class = d$class, id = NA_integer_, metric = "ca", value = 10,
               rep = rep, time = d$time, poly = "A"),
    data.frame(level = "class", class = d$class, id = NA_integer_, metric = "area_mn", value = 2,
               rep = rep, time = d$time, poly = "A")
  )
}

test_that("class areas a replicate-year lacks are written as 0 ha", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  aggRoot <- moduleFn(".nrvAggRoot", mod)
  fillClassAreaZeros <- moduleFn(".fillClassAreaZeros", mod)
  simRoot <- aggRoot(sim, "pm_zeros")
  ccRoot <- aggRoot(sim, "pm_zeros_CC")
  unlink(c(simRoot, ccRoot), recursive = TRUE)
  ## rep01 holds Decid in both years; rep02 Decid in year 700 and Pine in 800; today, Fir
  nrvtools::write_nrv_parquet(new_ca_replicate("rep01", "Decid", c(700L, 800L)), simRoot, "rep01")
  nrvtools::write_nrv_parquet(
    rbind(new_ca_replicate("rep02", "Decid", 700L), new_ca_replicate("rep02", "Pine", 800L)),
    simRoot,
    "rep02"
  )
  nrvtools::write_nrv_parquet(new_ca_replicate("rep01", "Fir", 0L), ccRoot, "CC")

  files <- fillClassAreaZeros(c(simRoot, ccRoot))

  v <- nrvtools::nrv_replicate_values(simRoot)
  ca <- v[v$metric == "ca", ]
  ## every replicate-year has the area of each class the subregion shows anywhere
  expect_equal(nrow(ca), 2L * 2L * 3L)
  expect_setequal(unique(ca$class), c("Decid", "Fir", "Pine"))
  expect_equal(ca$value[ca$replicate == "rep02" & ca$time == 800 & ca$class == "Decid"], 0)
  expect_equal(ca$value[ca$replicate == "rep01" & ca$time == 700 & ca$class == "Decid"], 10)
  expect_equal(sum(ca$value == 0), 12L - 4L)
  ## the other metrics stay undefined for an absent class
  expect_equal(nrow(v[v$metric == "area_mn", ]), 4L)
  ## today's map lacks Decid and Pine, which the replicates show
  cc <- nrvtools::nrv_replicate_values(ccRoot)
  expect_equal(cc$class[cc$metric == "ca"], c("Decid", "Fir", "Pine"))
  expect_equal(cc$value[cc$metric == "ca"], c(0, 10, 0))
  expect_setequal(basename(dirname(files)), c("replicate=rep01", "replicate=rep02", "replicate=CC"))

  ## a second call rebuilds the same zeros, not a second copy
  before <- lapply(files, function(f) as.data.frame(arrow::read_parquet(f)))
  expect_setequal(fillClassAreaZeros(c(simRoot, ccRoot)), files)
  expect_equal(lapply(files, function(f) as.data.frame(arrow::read_parquet(f))), before)
  expect_equal(nrow(nrvtools::nrv_replicate_values(simRoot)), nrow(v))
})

test_that("a patch metric's envelope spans the replicates' means over their patches", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  mod$ppFiles <- character(0)
  sim <- SpaDES.core::simInit(paths = testPaths)
  aggRoot <- moduleFn(".nrvAggRoot", mod)
  summarizeRepDataset <- moduleFn(".summarizeRepDataset", mod)
  root <- aggRoot(sim, "pm_perrep")
  unlink(root, recursive = TRUE)
  ## replicate r: patches of 1, 2 and 9 ha times r, and one patch whose age is undefined
  for (r in 1:3) {
    nrvtools::write_nrv_parquet(
      data.frame(level = "patch", class = "Pine", id = 1:4, metric = c(rep("area", 3), "sam_mdn"),
                 value = c(c(1, 2, 9) * r, NA), rep = sprintf("rep%02d", r), time = 700L,
                 poly = "A"),
      root,
      sprintf("rep%02d", r)
    )
  }

  env <- summarizeRepDataset(root, perReplicate = TRUE)
  area <- env[env$metric == "area", ]
  expect_equal(area$n_reps, 3L) ## replicates, not patches
  expect_equal(area$mean, mean(4 * 1:3))
  expect_equal(c(area$min, area$median, area$max), c(4, 8, 12))
  expect_false("sam_mdn" %in% env$metric) ## no replicate has a value
  expect_length(mod$ppFiles, 3L)

  pooled <- summarizeRepDataset(root)
  expect_equal(pooled$n_reps[pooled$metric == "area"], 9L) ## every patch of every replicate
})
