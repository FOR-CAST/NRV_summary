## The lw aggregates behind the leading-vegetation boxplots and large-patch histograms (2.0.0.9033):
## a vegetation type absent from a replicate-year gets 0, and each subregion gets a figure of the
## leading types' shares of its forest.

## leading proportions, shares of the forest and large-patch counts of `species` (and All species)
## in subregion `poly`, one replicate
new_lw_replicate <- function(rep, species, time, poly = "A") {
  d <- expand.grid(metric.1 = c(species, "All species"), class = c("Young", "Old"), time = time,
                   stringsAsFactors = FALSE)
  one <- function(level, metric, value) {
    data.frame(level = level, class = d$class, id = NA_integer_, metric = metric, value = value,
               rep = rep, time = d$time, poly = poly, metric.1 = d$metric.1)
  }
  rbind(one("landscape", "leadingProp", 0.5), one("landscape", "leadingShare", 0.1),
        one("class", "Npatch_ge100ha", 2))
}

test_that("vegetation types a replicate-year lacks get 0 in every lw metric", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  aggRoot <- moduleFn(".nrvAggRoot", mod)
  fillLeadingZeros <- moduleFn(".fillLeadingZeros", mod)
  simRoot <- aggRoot(sim, "lw_zeros")
  ccRoot <- aggRoot(sim, "lw_zeros_CC")
  unlink(c(simRoot, ccRoot), recursive = TRUE)
  ## subregion A: rep01 holds Decid in both years; rep02 Decid in 700 and Pine in 800; today, Fir.
  ## Subregion B only ever holds Pine.
  nrvtools::write_nrv_parquet(
    rbind(new_lw_replicate("rep01", "Decid", c(700L, 800L)),
          new_lw_replicate("rep01", "Pine", c(700L, 800L), poly = "B")),
    simRoot,
    "rep01"
  )
  nrvtools::write_nrv_parquet(
    rbind(new_lw_replicate("rep02", "Decid", 700L), new_lw_replicate("rep02", "Pine", 800L),
          new_lw_replicate("rep02", "Pine", c(700L, 800L), poly = "B")),
    simRoot,
    "rep02"
  )
  nrvtools::write_nrv_parquet(new_lw_replicate("rep01", "Fir", 0L), ccRoot, "CC")

  files <- fillLeadingZeros(c(simRoot, ccRoot))

  v <- as.data.frame(dplyr::collect(nrvtools::open_nrv_dataset(simRoot)))
  a <- v[v$poly == "A", ]
  ## each replicate-year has every metric and age class of each type subregion A shows anywhere
  expect_setequal(unique(a$metric.1), c("Decid", "Fir", "Pine", "All species"))
  expect_equal(nrow(a), 4L * 4L * 3L * 2L)
  expect_equal(nrow(unique(a[c("replicate", "time", "metric.1", "metric", "class")])), nrow(a))
  pick <- function(d, r, t, sp, met) d$value[d$replicate == r & d$time == t & d$metric.1 == sp &
                                               d$metric == met & d$class == "Old"]
  expect_equal(pick(a, "rep02", 800L, "Decid", "leadingProp"), 0)
  expect_equal(pick(a, "rep02", 800L, "Decid", "leadingShare"), 0)
  expect_equal(pick(a, "rep02", 800L, "Decid", "Npatch_ge100ha"), 0)
  expect_equal(pick(a, "rep01", 700L, "Decid", "leadingProp"), 0.5)
  expect_equal(unique(a$level[a$metric == "Npatch_ge100ha"]), "class")
  ## B never shows Decid or Fir, so it gets no rows for them
  expect_setequal(unique(v$metric.1[v$poly == "B"]), c("Pine", "All species"))
  ## today's map lacks Decid and Pine in A, which the replicates show
  cc <- as.data.frame(dplyr::collect(nrvtools::open_nrv_dataset(ccRoot)))
  cc <- cc[cc$poly == "A", ]
  expect_equal(pick(cc, "CC", 0L, "Pine", "leadingShare"), 0)
  expect_equal(pick(cc, "CC", 0L, "Fir", "leadingShare"), 0.1)
  expect_setequal(basename(files), "lw_zeros.parquet")

  ## a second call rebuilds the same zeros, not a second copy
  md5 <- unname(tools::md5sum(files))
  expect_setequal(fillLeadingZeros(c(simRoot, ccRoot)), files)
  expect_identical(unname(tools::md5sum(files)), md5)
  expect_equal(nrow(as.data.frame(dplyr::collect(nrvtools::open_nrv_dataset(simRoot)))), nrow(v))
})

## per-replicate shares of the forest: subregion A holds Pine and Bl_Spruce; B has none of either
new_shares <- function() {
  d <- expand.grid(rep = 1:3, class = c("Young", "Old"), metric.1 = c("Pine", "Bl_Spruce"),
                   poly = c("A", "B"), stringsAsFactors = FALSE)
  d$replicate <- sprintf("rep%02d", d$rep)
  d$time <- 700L
  d$metric <- "leadingShare"
  d$value <- ifelse(d$poly == "A", ifelse(d$metric.1 == "Pine", 0.3, 0.05), 0)
  rbind(d, transform(d, metric = "leadingProp", value = 0.5))
}

test_that("each subregion with leading forest gets one share-of-forest figure", {
  leadingShareFigTasks <- moduleFn(".leadingShareFigTasks", new.env())
  vegTypeLabels <- moduleFn(".vegTypeLabels", new.env())
  cc <- data.frame(poly = c("A", "B", "A"), class = "Old", metric.1 = "Pine",
                   metric = c("leadingShare", "leadingShare", "leadingProp"), value = 0.2)
  dir <- withr::local_tempdir()
  tasks <- leadingShareFigTasks(new_shares(), cc, c("Young", "Old"), dir, "WAU — ",
                                areaHa = function(poly) c(A = 1234.4, B = 10)[[poly]])

  expect_length(tasks, 1L) ## B has no leading forest of any type
  tk <- tasks[[1]]
  expect_equal(tk$plotter, "leadingShare")
  expect_equal(tk$file, file.path(dir, "A share of forest.png"))
  expect_equal(unique(tk$df$poly), "A")
  expect_equal(unique(tk$df$metric), "leadingShare")
  expect_equal(nrow(tk$cc), 1L)
  expect_equal(tk$title, "WAU — A")
  expect_match(tk$caption, "Total forested area in A: 1,234 ha", fixed = TRUE)
  expect_equal(tk$species_labels, stats::setNames(vegTypeLabels(c("Pine", "Bl_Spruce")),
                                                  c("Pine", "Bl_Spruce")))
})

test_that("a share-of-forest figure without a forested area has no area in its caption", {
  leadingShareFigTasks <- moduleFn(".leadingShareFigTasks", new.env())
  tasks <- leadingShareFigTasks(new_shares(), NULL, c("Young", "Old"), withr::local_tempdir(), "",
                                areaHa = function(poly) NA_real_)
  expect_no_match(tasks[[1]]$caption, "forested area")
  expect_null(tasks[[1]]$cc)
})

test_that("the share-of-forest figure renders to PNG", {
  mod <- new.env()
  leadingShareFigTasks <- moduleFn(".leadingShareFigTasks", mod)
  renderTasks <- moduleFn(".renderTasks", mod)
  ## render in this process: the worker count comes from P(sim), which needs the module loaded
  environment(renderTasks)$.plotWorkers <- function(sim, nTasks) 1L
  tasks <- leadingShareFigTasks(new_shares(), NULL, c("Young", "Old"), withr::local_tempdir(), "",
                                areaHa = function(poly) 100)
  files <- renderTasks(NULL, tasks)
  expect_equal(files, tasks[[1]]$file)
  expect_equal(unname(file.exists(files)), TRUE)
})
