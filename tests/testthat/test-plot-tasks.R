## The figure tasks plotFun() renders. The envelope figures carry the current condition
## (`<refCode>_CC`), and lm and pm add per-metric pages that set a layer's subregions side by side:
## both were dropped in the arrow rewire (2.0.0) and restored in 2.0.0.9031.

## an envelope of `metrics` x `polys` x two summary years, or its year-0 current-condition twin
new_envelope <- function(polys, metrics = c("ai", "ed"), cc = FALSE) {
  times <- if (cc) 0 else c(700, 800)
  d <- expand.grid(time = times, poly = polys, metric = metrics, stringsAsFactors = FALSE)
  d$class <- NA
  d$level <- "landscape"
  d$n_reps <- if (cc) 1L else 15L
  d$mean <- d$median <- d$min <- d$q25 <- d$q75 <- d$max <- seq_len(nrow(d))
  d
}

## one replicate of raw patch metrics: a patch-level metric (two patches per class and subregion)
## and a class-level one; `areaMn` adds the class's mean patch area
new_pm_replicate <- function(r, time = 700L, areaMn = FALSE) {
  patch <- expand.grid(id = 1:2, class = c("Decid", "Pine"), poly = c("A", "B"),
                       stringsAsFactors = FALSE)
  rbind(
    data.frame(level = "patch", class = patch$class, id = patch$id, metric = "area",
               value = patch$id * r, rep = r, time = time, poly = patch$poly),
    data.frame(level = "class", class = c("Decid", "Pine"), id = NA, metric = "ca",
               value = c(10, 20) * r, rep = r, time = time, poly = "A"),
    if (areaMn) {
      data.frame(level = "class", class = c("Decid", "Pine"), id = NA, metric = "area_mn",
                 value = 1.5 * r, rep = r, time = time, poly = "A")
    }
  )
}

## write the replicate and current-condition aggregates of one pm layer where plotFun reads them,
## and summarise them into `mod` as the pm event does
write_pm_aggregates <- function(sim, layer, mod, ...) {
  aggRoot <- moduleFn(".nrvAggRoot", mod)
  summarizeRepDataset <- moduleFn(".summarizeRepDataset", mod)
  refCode <- LandWebUtils::refCodeFor("pm", layer)
  unlink(c(aggRoot(sim, refCode), aggRoot(sim, paste0(refCode, "_CC"))), recursive = TRUE)
  for (r in 1:3) {
    nrvtools::write_nrv_parquet(new_pm_replicate(r, ...), aggRoot(sim, refCode), sprintf("rep%02d", r))
  }
  nrvtools::write_nrv_parquet(new_pm_replicate(1, time = 0L, ...),
                              aggRoot(sim, paste0(refCode, "_CC")), "CC")
  mod[[refCode]] <- summarizeRepDataset(aggRoot(sim, refCode), perReplicate = TRUE)
  mod[[paste0(refCode, "_CC")]] <- summarizeRepDataset(aggRoot(sim, paste0(refCode, "_CC")),
                                                       perReplicate = TRUE)
  invisible(refCode)
}

taskFiles <- function(tasks) vapply(tasks, function(tk) tk[["file"]], character(1))

test_that("each subregion's envelope figure carries that subregion's current condition", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  polys <- c("Boreal", "Foothills", "Montane")
  refCode <- LandWebUtils::refCodeFor("lm", "ANSR")
  mod[[refCode]] <- new_envelope(polys)
  mod[[paste0(refCode, "_CC")]] <- new_envelope(polys, cc = TRUE)
  envelopeFigTasks <- moduleFn(".envelopeFigTasks", mod)
  sim <- SpaDES.core::simInit(paths = testPaths)

  tasks <- envelopeFigTasks(sim, "lm", "ANSR", ylab = "value")

  expect_length(tasks, 12L) ## 2 metrics x 3 subregions x (ribbon, boxplot)
  expect_contains(basename(taskFiles(tasks)), c("Boreal ai_ribbon.png", "Montane ed_boxplot.png"))
  for (tk in tasks) {
    expect_equal(tk$cc$poly, unique(tk$df$poly))
    expect_equal(tk$cc$metric, unique(tk$df$metric))
  }
})

test_that("the per-metric pages set the subregions side by side, 12 to a page", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  polys <- sprintf("unit%02d", 1:13)
  refCode <- LandWebUtils::refCodeFor("lm", "ANSR")
  mod[[refCode]] <- new_envelope(polys)
  mod[[paste0(refCode, "_CC")]] <- new_envelope(polys, cc = TRUE)
  envelopeFigTasks <- moduleFn(".envelopeFigTasks", mod)
  sim <- SpaDES.core::simInit(paths = testPaths)

  tasks <- envelopeFigTasks(sim, "lm", "ANSR", ylab = "value", layout = "metric",
                            subdir = "per-metric", titlePrefix = "ANSR — ", minSubregions = 2L)

  files <- taskFiles(tasks)
  expect_setequal(
    basename(files),
    paste0(rep(c("ai", "ed"), each = 4), "_", rep(c("ribbon", "boxplot"), each = 2), "_p", 1:2, ".png")
  )
  expect_equal(unique(basename(dirname(files))), "per-metric")
  expect_equal(tasks[[1]]$title, "ANSR — ai (page 1 of 2)")
  expect_equal(tasks[[1]]$page, 1L)
  expect_setequal(tasks[[1]]$cc$poly, polys) ## every subregion's value, matched per panel
})

test_that("a layer with one subregion gets no per-metric pages", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  refCode <- LandWebUtils::refCodeFor("lm", "ANSR")
  mod[[refCode]] <- new_envelope("Boreal")
  envelopeFigTasks <- moduleFn(".envelopeFigTasks", mod)
  sim <- SpaDES.core::simInit(paths = testPaths)

  expect_length(
    envelopeFigTasks(sim, "lm", "ANSR", ylab = "value", layout = "metric", minSubregions = 2L),
    0L
  )
})

test_that("the per-class pages plot each replicate's values from the aggregates", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  write_pm_aggregates(sim, "FMA", mod)
  byClassFigTasks <- moduleFn(".byClassFigTasks", mod)

  tasks <- byClassFigTasks(sim, "pm", "FMA", xlab = "vegetation type")

  expect_setequal(
    basename(taskFiles(tasks)),
    c("area_box_p1.png", "area_violin_p1.png", "ca_box_p1.png", "ca_violin_p1.png")
  )
  area <- tasks[[match("area_box_p1.png", basename(taskFiles(tasks)))]]
  expect_equal(area$ylab, "Patch area (ha), mean over patches")
  ## one value per replicate, class and subregion: the mean of that replicate's two patches
  expect_equal(nrow(area$df), 3L * 2L * 2L)
  expect_setequal(area$df$value, c(1.5, 3, 4.5))
  expect_equal(nrow(area$cc), 4L) ## 2 classes x 2 subregions
  ca <- tasks[[match("ca_box_p1.png", basename(taskFiles(tasks)))]]
  expect_equal(ca$ylab, "Total class area (ha)")
})

test_that("a layer this run did not summarise gets no per-class pages", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  refCode <- write_pm_aggregates(sim, "Stale", mod)
  rm(list = refCode, envir = mod) ## its aggregates are on disk, but this run did not compute them
  byClassFigTasks <- moduleFn(".byClassFigTasks", mod)

  expect_length(byClassFigTasks(sim, "pm", "Stale"), 0L)
})

test_that("patch area has no per-class pages of its own beside area_mn", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  write_pm_aggregates(sim, "Pembina", mod, areaMn = TRUE)
  byClassFigTasks <- moduleFn(".byClassFigTasks", mod)

  files <- basename(taskFiles(byClassFigTasks(sim, "pm", "Pembina")))
  expect_setequal(files, paste0(rep(c("area_mn", "ca"), each = 2), "_", c("box", "violin"), "_p1.png"))
})

test_that("only units with a value get a panel, and the pages count those", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  polys <- sprintf("unit%02d", 1:13)
  refCode <- LandWebUtils::refCodeFor("lm", "ANSR")
  env <- new_envelope(polys)
  env[env$poly == "unit13", c("mean", "median", "min", "q25", "q75", "max")] <- NA
  mod[[refCode]] <- env
  envelopeFigTasks <- moduleFn(".envelopeFigTasks", mod)

  tasks <- envelopeFigTasks(sim, "lm", "ANSR", layout = "metric", subdir = "per-metric",
                            minSubregions = 2L)
  expect_setequal(
    basename(taskFiles(tasks)),
    paste0(rep(c("ai", "ed"), each = 2), "_", c("ribbon", "boxplot"), "_p1.png")
  )
  expect_false("unit13" %in% tasks[[1]]$df$poly)
  ## 12 panels: three rows of four, on a canvas of three rows
  expect_equal(c(tasks[[1]]$ncol, tasks[[1]]$nrow), c(4L, 3L))
  expect_equal(c(tasks[[1]]$width, tasks[[1]]$height), c(16, 2 + 2.75 * 3))
  ## seven panels take two rows
  mod[[refCode]] <- new_envelope(polys[1:7])
  tasks <- envelopeFigTasks(sim, "lm", "ANSR", layout = "metric", minSubregions = 2L)
  expect_equal(c(tasks[[1]]$nrow, tasks[[1]]$height), c(2, 2 + 2.75 * 2))
})

test_that("the figures name the metric and say what they show", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  polys <- c("Boreal", "Foothills")
  refCode <- LandWebUtils::refCodeFor("lm", "ANSR")
  mod[[refCode]] <- new_envelope(polys)
  mod[[paste0(refCode, "_CC")]] <- new_envelope(polys, cc = TRUE)
  envelopeFigTasks <- moduleFn(".envelopeFigTasks", mod)

  tasks <- envelopeFigTasks(sim, "lm", "ANSR")
  ribbon <- tasks[[1]]
  expect_equal(basename(ribbon$file), "Boreal ai_ribbon.png")
  expect_equal(ribbon$title, "Boreal — ai")
  expect_equal(ribbon$ylab, "Aggregation index (%)")
  expect_equal(ribbon$xlab, "Simulation year")
  expect_equal(
    gsub("\n", " ", ribbon$caption), ## wrapped to the canvas width
    paste(
      "Grey band: range of 15 replicates at each simulation year; line: mean; red dashed line:",
      "current condition (simulation year 0)."
    )
  )
  expect_equal(c(ribbon$width, ribbon$height), c(8, 5)) ## a single panel
  expect_match(
    gsub("\n", " ", tasks[[2]]$caption),
    "^Boxes: median and quartiles of 15 replicates at each simulation year; whiskers: minimum to maximum"
  )

  write_pm_aggregates(sim, "Caribou", mod)
  byClassFigTasks <- moduleFn(".byClassFigTasks", mod)
  tasks <- byClassFigTasks(sim, "pm", "Caribou", classLabels = function(x) paste0(x, "!"))
  box <- tasks[[match("area_box_p1.png", basename(taskFiles(tasks)))]]
  expect_setequal(box$df$class, c("Decid!", "Pine!"))
  expect_setequal(box$cc$class, c("Decid!", "Pine!"))
  expect_match(
    gsub("\n", " ", box$caption),
    paste0(
      "^Boxes: median and quartiles of 3 replicates x 1 simulation year; whiskers: .*",
      "red points: current condition \\(simulation year 0\\)\\. ",
      "A replicate's value is the mean over its patches\\.$"
    )
  )
  violin <- tasks[[match("ca_violin_p1.png", basename(taskFiles(tasks)))]]
  expect_match(violin$caption, "^Violins: distribution of 3 replicates x 1 simulation year")
})

test_that("the display names of vegetation types are LandWeb's group labels", {
  vegTypeLabels <- moduleFn(".vegTypeLabels", new.env())
  expect_equal(
    vegTypeLabels(c("Bl_Spruce", "Wh_Spruce", "Doug_fir", "Decid", "Mixed", "Pice_mar", NA)),
    c("Bl Spruce", "Wh Spruce", "Doug fir", "Decid", "Mixed", "Pice mar", NA)
  )
})

test_that("the figure tasks render to PNG", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  polys <- c("Boreal", "Foothills", "Montane")
  refCode <- LandWebUtils::refCodeFor("lm", "Caribou")
  mod[[refCode]] <- new_envelope(polys)
  mod[[paste0(refCode, "_CC")]] <- new_envelope(polys, cc = TRUE)
  write_pm_aggregates(sim, "Caribou", mod)
  envelopeFigTasks <- moduleFn(".envelopeFigTasks", mod)
  byClassFigTasks <- moduleFn(".byClassFigTasks", mod)
  renderTasks <- moduleFn(".renderTasks", mod)
  ## render in this process: the worker count comes from P(sim), which needs the module loaded
  environment(renderTasks)$.plotWorkers <- function(sim, nTasks) 1L

  tasks <- c(
    envelopeFigTasks(sim, "lm", "Caribou", ylab = "value")[1:2],
    envelopeFigTasks(sim, "lm", "Caribou", ylab = "value", layout = "metric",
                     subdir = "per-metric", minSubregions = 2L)[1:2],
    byClassFigTasks(sim, "pm", "Caribou")
  )
  files <- renderTasks(sim, tasks)

  expect_length(files, 8L)
  expect_setequal(files, taskFiles(tasks))
  expect_equal(unname(file.exists(files)), rep(TRUE, 8L))

  ## a figure that fails warns with its file name, and the others are still drawn
  bad <- tasks[[1]]
  bad$df <- bad$df[setdiff(names(bad$df), "max")]
  bad$file <- file.path(dirname(bad$file), "broken.png")
  expect_warning(files <- renderTasks(sim, list(bad, tasks[[2]])), "could not render .*broken[.]png")
  expect_equal(files, tasks[[2]]$file)
})
