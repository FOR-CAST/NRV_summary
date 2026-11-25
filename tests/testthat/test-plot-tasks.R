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
## and summarise them into `mod` as the pm event does. `times` are the years of each replicate
## (default: the one year of `...`); `codes = TRUE` writes the classes as the integer codes of an
## older aggregate, with their labels in a year-0 map at `mod$fvtm0`, and labels the envelopes
## as the pm event does. `unnumbered = TRUE` names the replicate folders "alpha", "beta", "gamma",
## so that the replicate number nrvtools parses from them (`rep`) is NA.
write_pm_aggregates <- function(sim, layer, mod, times = NULL, codes = FALSE, unnumbered = FALSE, ...) {
  aggRoot <- moduleFn(".nrvAggRoot", mod)
  summarizeRepDataset <- moduleFn(".summarizeRepDataset", mod)
  refCode <- LandWebUtils::refCodeFor("pm", layer)
  unlink(c(aggRoot(sim, refCode), aggRoot(sim, paste0(refCode, "_CC"))), recursive = TRUE)
  coded <- function(d) {
    if (codes) d$class <- c(Decid = 1L, Pine = 2L)[d$class]
    d
  }
  for (r in 1:3) {
    d <- if (is.null(times)) {
      new_pm_replicate(r, ...)
    } else {
      do.call(rbind, lapply(times, function(t) new_pm_replicate(r, time = t, ...)))
    }
    if (unnumbered) d$rep <- NA_integer_
    nrvtools::write_nrv_parquet(coded(d), aggRoot(sim, refCode),
                                if (unnumbered) c("alpha", "beta", "gamma")[r] else sprintf("rep%02d", r))
  }
  nrvtools::write_nrv_parquet(coded(new_pm_replicate(1, time = 0L, ...)),
                              aggRoot(sim, paste0(refCode, "_CC")), "CC")
  mod[[refCode]] <- summarizeRepDataset(aggRoot(sim, refCode), perReplicate = TRUE)
  mod[[paste0(refCode, "_CC")]] <- summarizeRepDataset(aggRoot(sim, paste0(refCode, "_CC")),
                                                       perReplicate = TRUE)
  if (codes) {
    vtm <- terra::rast(matrix(c(1, 2, 1, 2), 2), crs = "EPSG:3400")
    levels(vtm) <- data.frame(ID = 1:2, code = c("Decid", "Pine"))
    mod$fvtm0 <- file.path(tempfile("vtm"), "vegTypeMap_year0000.tif")
    dir.create(dirname(mod$fvtm0))
    terra::writeRaster(vtm, mod$fvtm0, datatype = "INT1U")
    for (rc in c(refCode, paste0(refCode, "_CC"))) {
      mod[[rc]] <- nrvtools::label_vegtype_classes(mod[[rc]], terra::rast(mod$fvtm0))
    }
  }
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

test_that("a run empties the figure folders of the events it draws, and only those", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  ppRoot <- moduleFn(".ppRoot", mod)
  clearFigDirs <- moduleFn(".clearFigDirs", mod)
  figs <- file.path(ppRoot(sim), "figures")
  withr::defer(unlink(figs, recursive = TRUE))
  old <- file.path(figs, c("boxplots/ANSR/Alpine Pice_gla.png", "histograms/BEC/100/Alpine Fir.png",
                           "lm/ANSR/Alpine ai_ribbon.png", "sspm/ANSR/Alpine seral.png"))
  for (f in old) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    file.create(f)
  }

  expect_setequal(clearFigDirs(sim, c("lw", "lm")), c("boxplots", "histograms", "lm"))
  expect_equal(unname(file.exists(old)), c(FALSE, FALSE, FALSE, TRUE))
  expect_length(list.files(file.path(figs, c("boxplots", "histograms", "lm")), recursive = TRUE,
                           include.dirs = TRUE), 0L)
})

## ---- envelopeStyle: a ribbon (the default) or one line per replicate and the median ------------------

## the figure of an envelope task, drawn as .renderTasks() draws it
plotOfTask <- function(tk) {
  nrvtools::plot_nrv_envelope(
    tk$df, type = tk$type, reps_df = tk$reps, facet = tk$facet, ylab = tk$ylab, title = tk$title,
    ncol = tk$ncol, nrow = tk$nrow, page = tk$page, cc = tk$cc, xlab = tk$xlab, caption = tk$caption
  )
}

## run plotFun() on `mod` for one layer, with `.renderTasks()` returning the tasks it was given
plotRun <- function(sim, mod, style, events = "pm", layer = "FMA") {
  mod$rptPolyNames <- layer
  plot <- moduleFn("plotFun", mod)
  env <- environment(plot)
  ## the helpers that read `mod` see this one
  for (f in c(".envelopeFigTasks", ".byClassFigTasks")) env[[f]] <- moduleFn(f, mod)
  env$P <- function(sim) {
    list(postprocessEvents = events, nrvWindow = 1, .studyAreaName = "", envelopeStyle = style)
  }
  env$.nrvPeriod <- function(sim) c(700, 800)
  got <- NULL
  env$.renderTasks <- function(sim, tasks) {
    got <<- tasks
    character(0)
  }
  plot(sim)
  got
}

layerClasses <- function(gg) unname(vapply(gg$layers, function(l) class(l$geom)[1L], character(1L)))

test_that("envelopeStyle defaults to ribbon: figures unchanged", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  polys <- c("Boreal", "Foothills")
  refCode <- LandWebUtils::refCodeFor("lm", "ANSR")
  mod[[refCode]] <- new_envelope(polys)
  mod[[paste0(refCode, "_CC")]] <- new_envelope(polys, cc = TRUE)
  write_pm_aggregates(sim, "FMA", mod)
  figTasks <- moduleFn(".envelopeFigTasks", mod)
  lmTask <- figTasks(sim, "lm", "ANSR", ylab = "value")[[1]]
  pmTask <- figTasks(sim, "pm", "FMA", classLabels = moduleFn(".vegTypeLabels", mod))[[1]]

  ## the layers and their data, recorded from the figures before envelopeStyle existed
  expect_equal(basename(lmTask$file), "Boreal ai_ribbon.png")
  gg <- plotOfTask(lmTask)
  expect_equal(layerClasses(gg), c("GeomRibbon", "GeomLine", "GeomHline"))
  d <- lapply(1:3, function(i) ggplot2::layer_data(gg, i))
  expect_equal(d[[1]]$x, c(700, 800))
  expect_equal(d[[1]]$ymin, c(1, 2))
  expect_equal(d[[1]]$ymax, c(1, 2))
  expect_equal(d[[2]]$y, c(1, 2))
  expect_equal(d[[3]]$yintercept, 1)

  expect_equal(basename(pmTask$file), "A ca_ribbon.png")
  gg <- plotOfTask(pmTask)
  expect_equal(layerClasses(gg), c("GeomRibbon", "GeomLine", "GeomHline"))
  d <- lapply(1:3, function(i) ggplot2::layer_data(gg, i))
  expect_equal(as.integer(d[[1]]$PANEL), 1:2)
  expect_equal(d[[1]]$ymin, c(10, 20))
  expect_equal(d[[1]]$ymax, c(30, 60))
  expect_equal(d[[2]]$y, c(20, 40))
  expect_equal(d[[3]]$yintercept, c(10, 20))

  ## through plotFun, "ribbon" builds the same tasks as not asking: no replicate values attached
  tasks <- plotRun(sim, mod, "ribbon")
  expect_setequal(basename(taskFiles(tasks)[vapply(tasks, function(t) t$plotter == "envelope", logical(1))]),
                  basename(taskFiles(figTasks(sim, "pm", "FMA",
                                              classLabels = moduleFn(".vegTypeLabels", mod)))))
  for (tk in tasks) {
    expect_null(tk$reps)
    expect_false(identical(tk$type, "lines"))
  }
})

test_that("lines style draws one line per replicate and the median", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  write_pm_aggregates(sim, "FMA", mod, times = c(700L, 800L))

  tasks <- plotRun(sim, mod, "lines")
  env <- Filter(function(t) t$plotter == "envelope", tasks)
  lines <- Filter(function(t) identical(t$type, "lines"), env)
  ca <- lines[[match("A ca_ribbon.png", basename(taskFiles(lines)))]]
  expect_equal(sort(unique(ca$reps$rep)), sprintf("rep%02d", 1:3)) ## the partitions

  gg <- plotOfTask(ca)
  expect_equal(layerClasses(gg), c("GeomLine", "GeomLine", "GeomHline")) ## no ribbon
  reps <- ggplot2::layer_data(gg, 1)
  expect_equal(nrow(reps), 12L) ## 3 replicates x 2 years x 2 classes
  for (pn in split(reps, reps$PANEL)) {
    expect_length(unique(pn$group), 3L) ## one line per replicate in each panel
    expect_equal(nrow(pn), 6L)
  }
  med <- ggplot2::layer_data(gg, 2)
  expect_equal(nrow(med), 4L) ## the median, one line per panel
  expect_equal(sort(med$y), c(20, 20, 40, 40))
  ## the boxplot figures are as before
  expect_true(all(vapply(Filter(function(t) identical(t$type, "boxplot"), env),
                         function(t) is.null(t$reps), logical(1L))))
})

test_that("lines style labels integer class codes from older aggregates", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  write_pm_aggregates(sim, "FMA", mod, times = c(700L, 800L), codes = TRUE)

  tasks <- plotRun(sim, mod, "lines")
  ca <- Filter(function(t) identical(t$type, "lines") && grepl("A ca_", basename(t$file)), tasks)[[1]]
  expect_setequal(ca$reps$class, c("Decid", "Pine")) ## the codes 1 and 2, labelled as the envelope's
  msgs <- testthat::capture_messages(gg <- plotOfTask(ca))
  expect_false(any(grepl("dropped", msgs))) ## every replicate row has a panel
  expect_equal(nrow(ggplot2::layer_data(gg, 1)), 12L)
})

test_that("replicate values are read once per dataset", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  write_pm_aggregates(sim, "FMA", mod, times = c(700L, 800L))
  refCode <- LandWebUtils::refCodeFor("pm", "FMA")
  roots <- character(0)
  real <- nrvtools::nrv_replicate_values
  testthat::local_mocked_bindings(
    nrv_replicate_values = function(x, ...) {
      roots <<- c(roots, as.character(x))
      real(x, ...)
    },
    .package = "nrvtools"
  )

  tasks <- plotRun(sim, mod, "lines")

  expect_true(any(vapply(tasks, function(t) identical(t$type, "lines"), logical(1L))))
  aggRoot <- moduleFn(".nrvAggRoot", mod)
  ## the replicates once, shared by the envelope figures and the per-class pages; the current condition once
  expect_equal(sum(roots == aggRoot(sim, refCode)), 1L)
  expect_equal(sum(roots == aggRoot(sim, paste0(refCode, "_CC"))), 1L)
})

test_that("an aggregate without replicate values falls back to the ribbon with a message", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  write_pm_aggregates(sim, "FMA", mod, times = c(700L, 800L))
  refCode <- LandWebUtils::refCodeFor("pm", "FMA")
  unlink(moduleFn(".nrvAggRoot", mod)(sim, refCode), recursive = TRUE)

  msgs <- testthat::capture_messages(tasks <- plotRun(sim, mod, "lines"))

  env <- Filter(function(t) t$plotter == "envelope", tasks)
  expect_gt(length(env), 0L)
  expect_setequal(vapply(env, function(t) t$type, ""), c("ribbon", "boxplot"))
  expect_equal(sum(grepl("no replicate values", msgs)), 1L) ## once, not once per figure
})

test_that("envelopeStyle rejects other values", {
  check <- moduleFn(".checkEnvelopeStyle", new.env())
  expect_identical(check("ribbon"), "ribbon")
  expect_identical(check("lines"), "lines")
  for (bad in list("dots", NA_character_, c("ribbon", "lines"), 1, NULL, character(0))) {
    expect_error(check(bad), "envelopeStyle must be \"ribbon\" or \"lines\"")
  }
})

test_that("envelopeStyle is a parameter that defaults to ribbon", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(md$parameters$default[[match("envelopeStyle", md$parameters$paramName)]], "ribbon")
})

test_that("lines figures render to PNG, with the NRV shading", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  write_pm_aggregates(sim, "FMA", mod, times = c(700L, 800L))
  renderTasks <- moduleFn(".renderTasks", mod)
  environment(renderTasks)$.plotWorkers <- function(sim, nTasks) 1L
  tasks <- Filter(function(t) identical(t$type, "lines"), plotRun(sim, mod, "lines"))
  tasks <- lapply(tasks, function(t) {
    t$nrvPeriod <- c(750, 800)
    t
  })
  expect_gt(length(tasks), 0L)

  files <- renderTasks(sim, tasks)

  expect_setequal(files, taskFiles(tasks))
  expect_true(all(file.exists(files)))
  ## the shading, behind the lines: a rectangle over the NRV years, on the continuous time axis
  gg <- moduleFn(".addNrvShading", mod)(plotOfTask(tasks[[1]]), c(750, 800))
  rect <- Filter(function(l) inherits(l$geom, "GeomRect"), gg$layers)
  expect_length(rect, 1L)
  expect_identical(gg$layers[[1]], rect[[1]]) ## drawn first, behind the data
  d <- ggplot2::layer_data(gg, 1)
  expect_equal(unique(c(d$xmin, d$xmax)), 800) ## the NRV years that the figure has (year 800)
})

test_that("lines are grouped by replicate, not by the replicate number parsed from its folder", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  write_pm_aggregates(sim, "FMA", mod, times = c(700L, 800L), unnumbered = TRUE)

  tasks <- plotRun(sim, mod, "lines")
  ca <- Filter(function(t) identical(t$type, "lines") && grepl("A ca_", basename(t$file)), tasks)[[1]]
  reps <- ggplot2::layer_data(plotOfTask(ca), 1)

  for (pn in split(reps, reps$PANEL)) {
    expect_length(unique(pn$group), 3L) ## three lines, though no folder gives a number
  }
})

test_that("lm envelopes draw as lines too: the subregion figures and the per-metric pages", {
  withr::local_package("SpaDES.core")
  mod <- new.env()
  sim <- SpaDES.core::simInit(paths = testPaths)
  polys <- c("A", "B")
  refCode <- LandWebUtils::refCodeFor("lm", "FMA")
  root <- moduleFn(".nrvAggRoot", mod)(sim, refCode)
  unlink(root, recursive = TRUE)
  one <- function(r) {
    d <- expand.grid(time = c(700L, 800L), poly = polys, metric = c("ai", "ed"), stringsAsFactors = FALSE)
    data.frame(level = "landscape", class = NA_character_, id = NA, metric = d$metric, value = r * d$time / 100,
               rep = r, time = d$time, poly = d$poly)
  }
  for (r in 1:3) nrvtools::write_nrv_parquet(one(r), root, sprintf("rep%02d", r))
  mod[[refCode]] <- moduleFn(".summarizeRepDataset", mod)(root, perReplicate = TRUE)

  tasks <- plotRun(sim, mod, "lines", events = "lm")

  lines <- Filter(function(t) identical(t$type, "lines"), tasks)
  expect_length(lines, 2L * 2L + 2L) ## per subregion x metric, and one page per metric
  expect_setequal(basename(taskFiles(lines)),
                  c("A ai_ribbon.png", "B ai_ribbon.png", "A ed_ribbon.png", "B ed_ribbon.png",
                    "ai_ribbon_p1.png", "ed_ribbon_p1.png"))
  page <- lines[[match("ai_ribbon_p1.png", basename(taskFiles(lines)))]]
  reps <- ggplot2::layer_data(plotOfTask(page), 1)
  expect_equal(nrow(reps), 12L) ## 3 replicates x 2 years x 2 subregions
  for (pn in split(reps, reps$PANEL)) expect_length(unique(pn$group), 3L)
  expect_true(all(vapply(Filter(function(t) identical(t$type, "boxplot"), tasks),
                         function(t) is.null(t$reps), logical(1L)))) ## no slices built for the boxplots
})
