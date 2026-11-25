Known issues: <https://github.com/FOR-CAST/NRV_summary/issues>

# NRV_summary (development version)

* **On the `cc` figures, a unit named "NA" is no longer the missing unit (`2.0.0.9047`).** Rows were grouped into figures by `paste(layer, poly)`, which writes a missing unit and a unit named "NA" alike, so their rows shared one figure. The grouping key now keeps them apart; nothing is written differently otherwise.

* **`.autocorrSummary()` no longer guesses the number of snapshots per replicate (`2.0.0.9047`).** The default, from the length of the NRV period and `summaryInterval`, ignored that the summary grid starts at `summaryPeriod[1]`. `nSnapshotsPerRep` is now required; the `autocorr` event already passed the exact count, so its results are unchanged.

* **Tests for the stop on a missing year-0 cohort file (`2.0.0.9047`).** The `st` event and the `bc` year-0 check stop with both formats' paths named when the current condition's cohort file is missing; this had no test. No change to the module.

* **A cohort file is paired with its `pixelGroupMap` by replicate folder and year, however the `.qs2` year is padded (`2.0.0.9047`).** Mode 'multi' used to build the `.qs2` name from the `pixelGroupMap` name, so `cohortData_year700.qs2` beside `pixelGroupMap_year0700.tif` was not found and the `bc` and `st` events stopped. Now, among the cohort files `InitMulti()` found, the one of the same replicate folder and the same year (compared as a number) is used, by every site that chooses a cohort file. The one named as the caller names it is taken if it is among them, else the first by name, with a warning from `InitMulti()` when a folder and year have several. Parquet part files are still matched by their exact `year=<YYYY>` path. With both formats the parquet file still wins. Files padded as LandWeb's are (`cohortData_year0700.qs2`) are chosen as before.

* **`cohortData` can be saved as a parquet dataset grouped by simulation year (opt-in), and mode 'multi' reads both formats (`2.0.0.9046`).** Mode 'single' still saves each replicate's `cohortData` (at the start, at each summary year and at the end) as `cohortData_year<YYYY>.qs2` by default, and saves and registers the same files as before. With the new parameter `cohortDataFormat = "parquet"`, it saves them with `nrvtools::write_cohort_parquet()` instead, as `cohortData/year=<YYYY>/part-0.parquet` in the replicate folder (FOR-CAST/NRV_summary#1). `cohortDataFormat` is `"qs2"` (default) or `"parquet"`; any other value stops `init`, in either mode. The `.qs2` files are about a quarter of the size and faster to write and read; parquet is an open format that other software reads, and the years of a replicate can be queried together with `arrow::open_dataset()` (the module documentation gives measurements). Mode 'multi' reads both formats with `nrvtools::read_cohort_data()`, so runs saved before this version are summarised as before. A cohort file counts when it was found (listed in the `outputsDF`, or on disk without one) and exists. `InitMulti()` takes a replicate and year as complete when a file of either format counts (the message for a missing one names both paths) and uses the parquet file when both do, with a warning that gives the number of such replicate-years and names up to three of their parquet files (neither file is removed). Without an `outputsDF`, it lists the `cohortData/` folder of every folder in `outputPath` (no other folder is listed further down) and keeps the part files of the replicates. Each `pixelGroupMap` is paired with the cohort file of the same folder and year; when the `bc` or `st` event, the two that read these pairs, will run, a `pixelGroupMap` without one stops `InitMulti()`, naming it. The reporting-group maps (`sppEquivColReporting`), the `bc` seral maps and the `st` event read either format and choose by the same rule, so all read the same file of a replicate and year; a missing year-0 (current-condition) cohort file, which `InitMulti()` does not check, stops them with both paths named.

* **nrvtools 0.5.0 or later is required (`2.0.0.9046`).** `write_cohort_parquet()` and `read_cohort_data()` are in 0.5.0.

* **LandWeb's existing outputs are summarised as before (`2.0.0.9046`).** On InterlakePlain (15 replicates, `.qs2` cohorts, no `outputsDF`, LandWeb's settings), `InitMulti()` finds and selects the same files as `2.0.0.9045`, and the reporting-group maps it builds from the cohorts are byte-identical; `read_cohort_data()` returns the same table as `qs2::qs_read()` for each of the 120 cohort files.

* **`InitMulti()` drops the names of `dirnamesFromSet()`'s result (`2.0.0.9046`).** PredictiveEcology/SpaDES.core#543 names the replicate folders it returns by replicate id; the folders are the same, and `InitMulti()` now gives the same `allReps` and first replicate folder with either version.

* **`envelopeStyle = "lines"` draws the lm and pm envelopes as one line per replicate (`2.0.0.9045`).** The new parameter `envelopeStyle` is `"ribbon"` (default) or `"lines"`; any other value stops `init`. With `"lines"`, each lm and pm envelope figure (the per-subregion figures of both, and the lm per-metric pages; the pm per-metric pages are the per-class ones, which are unchanged) shows a thin grey line for each replicate under a bold line for their median, drawn by `nrvtools::plot_nrv_envelope(type = "lines")`; the file names stay `<subregion> <metric>_ribbon.png` and `<metric>_ribbon_p<N>.png`, and the caption says what the lines are. A patch-level metric's line is each replicate's mean over its patches, as in the per-class pages. The replicate values are read once per reporting layer and shared with the pm per-class pages, and their class codes are labelled as the envelope's are, so integer codes in an older aggregate match their panels. A layer whose aggregates hold no replicate values gets ribbons, and a message says so. The boxplot figures, the `bc` seral figures and the per-class pages are drawn the same in both styles, and the NRV shading still applies. With `"ribbon"` the figures are drawn as before.

* **Opt-in `cc` event: the lm and st metrics against the NRV, on one figure per unit (`2.0.0.9045`).** When `"cc"` is in `postprocessEvents` (not in the default; mode 'multi'), with `"lm"` or `"st"` in the same run, the event draws, for each reporting layer and unit, one figure with the current condition of the lm metrics and the stand-level (st) measures against their NRV ranges (the pm, lw and bc metrics are not on it), each range rescaled to 0-1: `figures/all/<layer>/<unit>.png`. Layer folders and unit files are named so that two names never share a file (`A B` and `A_B` get two) and none can leave `figures/all/`; the folder is emptied when the event runs and the `plot` event leaves it alone. The stand-level rows (from `st`) go on the figure of every unit named like the study area (`.studyAreaName`); with no such unit they have a figure of their own, `figures/all/<studyArea>/<studyArea>.png`. Age classes run from young to old, in the order of `ageClasses`. The event writes no table itself: with `"cc"` in the run, the lm event writes `csv/lm/<layer>/lm_nrvComparison.csv`. It runs after lm, st and autocorr and before the plot event.

* **Opt-in `st` event: stand-level measures against the NRV (`2.0.0.9045`).** When `"st"` is in `postprocessEvents` (not in the default; mode 'multi'), the module reads the replicates' saved snapshots at each summary time and takes the share of the forest in each age class, the forested share of the area, the share of the area led by each type (the mixed class included) and the total biomass. It writes `csv/st/<studyArea>/st.csv`, `st_CC.csv` (the current condition) and `st_nrvComparison.csv` (whether the current condition is inside, above or below the NRV range), and figures in `figures/st/<studyArea>/` (`stand_age.png`, `stand_leading.png` (drawn only when some type leads at least 0.5% of the area, in the replicates or today), `stand_biomass.png`, `species/lead_<code>.png` and the paged `current_condition_overview_p<N>.png`), whose folder is emptied when the event runs. With `"lm"` in the same run, and `"st"` or `"cc"` among the events, `csv/lm/<layer>/lm_nrvComparison.csv` is also written, from the landscape-level rows. With `"stability"` in the same run, the stand series are tested too. Differences from the PredictiveEcology fork: the biomass sums forested cells only (urban-masked cells excluded), the last page of the overview draws only its own panels, a code shared by several `sppEquiv` rows is named by all of them, and the comparison table keeps landscape-level rows only.

* **nrvtools 0.4.0 or later is required (`2.0.0.9045`).** `plot_nrv_envelope(type = "lines")` is in 0.4.0.

* **No default changes (`2.0.0.9045`).** With the defaults (`postprocessEvents = c("lm", "pm")`, `envelopeStyle = "ribbon"`), no event, file or figure is added or changed.

* **Annual landscape series and autocorrelation check, both opt-in (`2.0.0.9044`).** With `recordAnnualSeries = TRUE` (default `FALSE`; mode 'single'), each replicate writes `annualSeries.csv` for the reporting area: the proportion of forest young and old, the mean stand age, the total biomass and the proportion led by each species, one row per year (about 1,000 extra events per replicate). With `"autocorr"` in `postprocessEvents` (not in the default; mode 'multi'), the module takes each series' autocorrelation over the NRV years, around the mean of all replicates with the variance pooled over them, and writes `csv/autocorrelation/autocorrelation.csv` and `figures/autocorrelation/all/`: the autocorrelation time (that of the exponential decay with the lag-1 autocorrelation as its rate; no curve is fitted), the correlation between snapshots `summaryInterval` years apart, the effective number of independent snapshots and a recommended thinning interval. The estimates are lower bounds, because a short series biases the sample autocorrelation toward 0; the event says so when the NRV has fewer than 20 times the recommended thinning in years. If the replicates' means differ, the common mean errs toward more thinning. The longest lag is the largest of `autocorrMaxLag` (200 years), `summaryInterval` and 50, cut to the NRV years less one. The event says how many replicates it read, leaves out (with a warning that names the replicate) a file that does not cover the NRV years, and names the series it leaves out (missing from a replicate, constant, or with fewer than 3 values). Needs LandR 1.2.0.9054 or later (`leadingByPixelGroup()`). With the defaults, no event, file or figure is added.

* **`nrvWindow`: the NRV can be the last part of `summaryPeriod` (`2.0.0.9043`).** `nrvWindow` is the fraction of `summaryPeriod` that counts as the natural range of variation (default 1, the whole period). The lw pooling, the pm by-class pages and the NRV shading use those years, and the shading is drawn only when `nrvWindow` is below 1. Reused aggregates are recomputed when the NRV years change; as existing lw aggregates carry no stamp of their years, the first run with `reuseAggregates = TRUE` after this update rebuilds them once, with the same results. `init` stops with a plain message when no summary time falls in the NRV years.

* **Opt-in `stability` event (`2.0.0.9043`).** When `"stability"` is in `postprocessEvents`, the module tests each lm and pm metric series of the events run alongside it (and the burn series, read from burnSummaries' `burnSummaries_fireBlocks.csv`; a message says so if it is absent) for a trend, writes `csv/stability/all/stability.csv`, and draws paged figures with each page showing only its own series. New parameters: `stabilityAlpha`, `stabilityMinChange` and `stabilityWindow`. A message says so when the run has no `lm` or `pm` event; `init` warns when the stability window holds fewer than 3 summary times (e.g. the defaults, `summaryInterval = 100`), and the verdict then says that no series could be tested rather than that the run is long enough.

* **Current-condition overview (`2.0.0.9043`).** When `"stability"` and `"lm"` are both in `postprocessEvents`, the stability event also draws where each lm metric's current condition sits against its NRV range, on each metric's own scale, up to 24 panels per page.

* **LandWeb summaries are unchanged (`2.0.0.9043`).** LandWeb leaves `nrvWindow` at 1 and does not run `stability`. On WesternBoreal (15 replicates), all 192 lw CSVs are byte-identical to those of `2.0.0.9042`; the only new files are the `nrvYears.txt` stamps beside each layer's lw aggregates.

* **`vegTypeMap` is saved as INT1U, and a `timeSeriesTimes` year off the `summaryInterval` grid no longer starts a second chain (`2.0.0.9042`).** Every replicate's saved `vegTypeMap_year*.tif` is now 8-bit (it keeps its colour table; values are unchanged), and the run stops if a code of 255 or more would not fit (255 is the NA value of an 8-bit raster). A `timeSeriesTimes` year inside `summaryPeriod` but off the grid (e.g. 725 for `summaryPeriod = c(700, 1000)`, `summaryInterval = 50`) runs once and no longer schedules extra years (775, 825, ...).

* **`timeSeriesTimes = NA` takes the default (`2.0.0.9042`).** The default is `start + 601:650` when the run is that long, otherwise the last 50 years of the run. An explicit `NA` used to stop `init` with a missing-value error.

* **`mode = "multi"` can take the replicate files from an `outputsDF` input (`2.0.0.9041`).** `outputsDF` is a new optional input: a data.frame with a `file` column, such as the replicate simulations' `outputs(sim)` tables row-bound together. When it is supplied, `InitMulti()` reads those files, and the replicate folders (`rep1`, `rep2`, ..., which may lie outside `outputPath`) are found with `SpaDES.core::dirnamesFromSet()`. Without it, the replicates are the folders under `outputPath` as before (`rep01`, `rep02`, ...), and `rep01`-style folders now also work when an outputs table is given.

* **`InitMulti()` checks the replicate files and the replicate count (`2.0.0.9041`).** Missing files are reported in both cases (with or without `outputsDF`; a file listed in `outputsDF` must also exist on disk), and the run stops if fewer replicates are found than the `reps` parameter asks for. The cohort, pixelGroup and map file filters now match on the file name, not on the whole path.

* **`simTimes` defaults to the run's start and end (`2.0.0.9041`).** When `simTimes` is not set, `InitMulti()` takes it from the simulation's times (`SpaDES.core::resolveSimYears()`), so runs that start after year 0 find their year-0 rasters. A `simTimes` that is set is used as given.

* **BC seral stages and the stand-age animation read the first replicate found (`2.0.0.9041`).** They no longer assume a folder called `rep01`.

* **SpaDES.core 3.2.1.9057 or later is required (`2.0.0.9041`).** `InitMulti()` uses `dirnamesFromSet(leafOnly = FALSE)` and `resolveSimYears()`, relies on `dirnamesFromSet()` keeping zero-padded folders such as `rep01` (from 3.2.1.9038), and on its bare-number folder fallback (`1/`, `2/`; from 3.2.1.9057).

* **LandWeb multi-mode summaries are unchanged (`2.0.0.9041`).** For LandWeb (`rep01` folders, no `outputsDF`), `InitMulti()` selects the same files as before (checked on one study area's outputs).

* `NRV_summary.md` shows the module's events, the objects each reads and writes, and a figure per event, generated from the code (needs SpaDES.core 3.2.1.9068 or later; older versions leave them out).

* **A run shorter than 650 years gets a time series (`2.0.0.9040`).** The default `timeSeriesTimes` is the last 50 years of the run when the run ends before year 650, and years 601 to 650 otherwise.

* **A study area with no tree species is skipped (`2.0.0.9040`).** If `sppEquiv` holds no tree species, the module summarises and saves nothing and says so in a message.

* **Each summary year is summarised once (`2.0.0.9040`).** Year 1000 was summarised and saved twice; every summary time now gets one summary and one save.

* **Inputs from another module's init are accepted (`2.0.0.9040`).** In `mode = "single"` the inputs are checked at the first summary event rather than at init, so a module loaded after NRV_summary can supply them; a missing input is reported by name, all at once.

* **`InitMulti()` matches analysis years and replicate folders exactly (`2.0.0.9040`).** In `mode = "multi"`, a year is read from the `_year<YYYY>.` in a file name and a replicate from a whole folder name, so digits elsewhere in a path no longer select every year, and a folder elsewhere in the path that contains a replicate's name (e.g. `rep01_old`) or a run of 100 or more replicates (`rep10` matching `rep100`) no longer selects the wrong files. Auxiliary `.aux.xml` files are no longer read as cohort or pixel-group files. On LandWeb's WesternBoreal outputs (15 replicates, years 0 to 1000) it selects the same files as before, so its multi-mode summaries are unchanged.

* **The module declares every package it calls (`2.0.0.9040`).** `arrow`, `crayon`, `future`, `reproducible` and `rlang` are added to `reqdPkgs`, and `pemisc` requires 0.0.4.9016 or later.

* **A vegetation type absent from a snapshot has no leading proportion there (`2.0.0.9036`).** `leadingProp` is the share of a type's own forest in each age class, which has no value (0 of 0 ha) where the type is absent, so the `lw` zero rows (`lw_zeros.parquet`) no longer fill it: the leading-proportion boxplots and CSVs summarise the replicate-years in which the type leads. The area metrics (`leadingShare`, large-patch counts and areas) keep their zeros, and a type absent today has no current-condition `leadingProp`. **This changes the leading-proportion boxplots and CSVs** wherever a type is absent from some replicate-years (in LandWeb's WesternAlbertaUpland, fir in all 10 tenures).

* **A run removes the figures of earlier runs before drawing its own (`2.0.0.9035`).** In `mode = "multi"`, `figures/lm`, `pm`, `boxplots`, `histograms` and `sspm` are emptied before the events that draw into them; the folders of events not run keep their figures. Figures of a type, subregion or reporting layer that was renamed or dropped stayed beside the current ones: after the 10-22 WesternAlbertaUpland rebuild, 5,749 of 16,013 figures were from earlier runs, such as per-species boxplots named by species code. The current condition in the leading-type boxplots and the share-of-forest figures is a red circle, so simulated values under it stay visible (nrvtools >= 0.3.4.9003).

* **The current condition's stand ages are today's, not the binned ages of the cohorts (`2.0.0.9034`).** In `mode = "single"`, the year-0 `standAgeMap` now comes from the age layer `sim$standAgeMap` (Biomass_borealDataPrep's, cropped to `pixelGroupMap`), on the pixels with cohorts. It came from the year-0 cohorts, whose ages Biomass_borealDataPrep bins (on WesternAlbertaUpland, ages above 70 to the nearest 20 years): there, 30.4% of the forest in the reporting area was 120 years or older by the cohorts, against 22.1% by the age layer. Each age is rounded half up and then floored to the decade, as `LandR::standAgeMapGenerator()` floors the stand ages of every simulated year, so the current condition has the 10-year resolution of the envelopes it is compared with. The flooring leaves the age-class shares unchanged (the class limits are multiples of 10) and lowers the patch ages: on WesternAlbertaUpland, the current condition's patch median age (`sam_mdn`), averaged over the patches of the reporting area, is 71.8 years, against 74.5 by the cohorts and 76.3 by unfloored ages. A pixel with cohorts but no value in the age layer keeps its cohort age. The simulated years keep the cohorts' ages. Of the summaries, only the current condition's patch ages (`pm`, `sam_mdn`) read this map; the `lw` summaries use time since fire. New parameter `ccStandAgeSource`: `"standAgeMap"` (the default) or `"cohortData"`, the former map; without `sim$standAgeMap` on a grid aligned with `pixelGroupMap`, the module warns and uses the cohorts. **This changes the year-0 stand-age map of every replicate and the current condition's `pm` age values.**
* **A figure of the leading vegetation types' shares of each subregion's forest, and absent types count as 0 in the `lw` summaries (`2.0.0.9033`).** `boxplots/<layer>/<subregion> share of forest.png` has a panel per age class and a box per leading type of its share of the subregion's forested area (`leadingShare`, from nrvtools >= 0.3.4.9002), with the current condition in red. Every box has the same denominator, so a type that leads on a few pixels has a small share, where its proportions by age class swing between 0 and 1. The nrvtools producers write rows only for the types present in a map, so a type that led in 7 of 105 replicate-years had its boxes and histograms drawn from those 7 alone (Fir in WAU's Subalpine ANC: median Old proportion 1). The `lw` aggregates now hold 0, in every metric, for each type a subregion shows (in its current condition or in any replicate-year) but lacks in a replicate-year, as the `pm` class areas do. In the leading-proportion boxplots such a replicate-year counts as 0 in every age class, so there the type's age classes sum to 0 rather than 1.
* **Patch metric envelopes span the replicates, and class areas count absent classes as 0 ha (`2.0.0.9032`).** The `pm` envelopes of a patch-level metric (`area`, `sam_mdn`) pooled every patch of every replicate, so the band ran from a single cell to the largest patch, and the current condition's patch mean and median, drawn on the ribbon and box figures, were up to 51 times apart. Each replicate, and the current condition, now contributes its mean over its patches (`nrvtools::nrv_replicate_values()`), as on the per-class pages; the envelope CSVs change accordingly, and their `n_reps` counts replicates, not patches. For `lm` and `pm`, `n_reps` counts the replicates with a value: one whose metric is undefined there (the standard deviation of a single patch, say) is left out. The per-class pages of `area` are dropped: a replicate's mean patch area is the class's mean patch area, `area_mn`, which has its own. landscapemetrics reports no class area for a class absent from a map, so the `ca` boxes and envelopes summarised only the replicate-years where the class occurs. After the `pm` aggregation, a class area of 0 is written (as `replicate=<id>/ca_zeros.parquet`) for each class that a subregion shows in its current condition or in any replicate-year but lacks in a replicate-year, the current condition included. Patch size, shape and age metrics stay undefined for an absent class. **This changes the `pm` envelopes, CSVs and figures.**
* **The NRV figures say what they show (`2.0.0.9032`).** The y axis gives the metric's plain name and unit (`nrvtools::nrv_metric_info()`; a patch-level metric adds "mean over patches") with thousands separators, the x axis of the envelopes reads "Simulation year", vegetation types are shown by their display names ("Bl Spruce"), and each figure has a caption built from its data: on the ribbons, "Grey band: range of <n> replicates at each simulation year; line: mean; red dashed line: current condition (simulation year 0)"; on the envelope boxplots, the median and quartiles of the replicates with whiskers from minimum to maximum; on the per-class pages, the median and quartiles of <n> replicates x <k> simulation years with whiskers to 1.5 x IQR, or their violins. The titles are unchanged (subregion or layer, then metric). A page has as many rows of four panels as it needs (up to three), with the canvas height to match, a single-panel figure a smaller canvas, and every envelope and per-class figure the same strip text size. The envelope quartiles are exact (nrvtools >= 0.3.4.9001), where they were Arrow's approximations.
* **One figure that fails no longer stops the others (`2.0.0.9032`).** A figure that fails to render warns with its file name. The pages of a metric are counted from the subregions (and classes) that have a value, the panels drawn, where a layer of 13 or more subregions could ask for a page that did not exist. The per-class pages are drawn only for a layer this run summarised, rather than from whatever aggregates were on disk. Two runs over the same `lm` and `pm` aggregates write identical envelope CSVs, class-area zero files and figures: the zero rows are sorted, and the per-replicate means no longer depend on the order in which Arrow reads the rows. Requires nrvtools >= 0.3.4.9001.
* **The envelope figures show the current condition, and the per-metric pages are back (`2.0.0.9031`).** Every `lm`, `pm` and `bc` envelope figure draws the current condition (the year-0 snapshot, `<refCode>_CC`) as a red dashed line in each panel where it has a value, matched by subregion, class and metric: the current-condition mean on the ribbon figures, its median on the box-and-whisker figures. The `lm` and `pm` analyses again write, for each metric, figures that set a layer's subregions side by side, 12 panels to a page, under `figures/<kind>/<layer>/per-metric/` (`<metric>_<type>_p<page>.png`). For `lm` these are the envelopes over time, a panel per subregion, skipped for a layer with one subregion. For `pm` they are box and violin plots of each vegetation type's values, a panel per subregion, with the current condition as a red point; the values are those of each replicate and summary year, a patch-level metric averaged over each replicate's patches. A vegetation type present today but in no replicate or summary year, or whose metric is undefined in all of them, keeps its point, without a box. Both views were dropped in the arrow rewire (`2.0.0`), where the box and violin plots had shown the across-replicate means of the summary years. Requires nrvtools >= 0.3.4.9000.
* **Year 0, the current condition, is saved before any time-0 dynamics (`2.0.0.9030`).** It was saved at `.last()`, after Biomass_core's growth and mortality when its `growthInitialTime` is `start(sim)`, so today's values held a simulated year of growth and old-age mortality and moved with the species traits (on WesternAlbertaUpland, that year changed the leading type of 4.4% of pixels). A `timeSeriesTimes` or `summaryPeriod` value equal to `start(sim)` no longer adds a second, later year-0 save. **This changes the year-0 files of every replicate and so the current-condition values.**
* **New parameter `sppEquivColReporting` reports vegetation types by groups of the simulated species (`2.0.0.9029`).** In `mode = "multi"`, the maps that the `lm`, `lw` and `pm` analyses read, for the simulated years and for current conditions, are rebuilt from the saved `cohortData` and `pixelGroupMap` with each species recoded to its group in that column of `sppEquiv`, and written under `postprocess/_reportingMaps/`. `LandR::vegTypeMapGenerator()` sums biomass only within a species code, so this differs from relabelling the species-level maps: a stand of 40% black spruce, 40% tamarack and 20% aspen is led by the group holding both conifers. A simulated code in two groups, or a group holding both conifers and broadleaves, stops the run. Colours come from the new optional input `sppColorVectReporting`. `NA`, the default, keeps the species-level maps. The grouping is `LandWebUtils::sppEquiv_groups()` and `recode_cohorts()` (LandWebUtils >= 1.0.3.9045), shared with LandMine.
* **`vegLeadingProportion` defaults to `LandR::leadingSpeciesProp()` (`2.0.0.9028`),** as in Biomass_core: option `LandR.leadingSpeciesProp`, which takes `LandR.mixedwoodProp` (0.75) unless set. **The default changes from 0.8 to 0.75, which changes the vegetation type maps** of a run that does not set it. With the old default, `paramCheckOtherMods()` stopped a run whose Biomass_* modules used 0.75. Requires LandR >= 1.2.0.9024.
* **`postprocessEvents` is checked at init (`2.0.0.9027`).** `'fd'` was scheduled as a stub that stopped in `browser()` after the other summaries had run, and any other unknown value was skipped without a word. Both now stop the run before any work, naming the values at fault.
* **The metric events set their worker count without a warning or an error (`2.0.0.9027`).** They tweaked `workers` on the caller's plan inline: under `sequential` that warned, and under a plan set as `plan(<backend>, workers = n)` it stopped with 'formal argument "workers" matched by multiple actual arguments', because future (<= 1.76.0) appends a tweak instead of replacing it. The BC seral event called `plan(workers = n)`, which changed nothing. All four now use `.planWithWorkers()`.
* **The current-condition urban correction now completes (`2.0.0.9026`).** `.maskCC()` wrote its masked copies as `_cc/cc_vegTypeMap.tif` etc., but nrvtools labels each per-map result from the path and requires `<rep>/<prefix>_year<YYYY>.tif`, so the landscape-metric summaries stopped with "cannot parse metric label(s)" the first time the correction ran (it had always been skipped before, because `LandTypeCC_reporting` was missing or on the wrong grid). The copies now keep the source's `<rep>/<file>` tail under `_cc/`, so their labels match the unmasked maps'.
* **`mode = "multi"` outputs are now registered with `registerOutputs()` (`2.0.0.9025`).** None were, so a targets stage running this module tracked 0 files: its outputs could change or disappear without the pipeline noticing, while `burnSummaries` tracked 557 MB for the same run. The parquet aggregates, envelope CSVs and stand-age animation are recorded as they are written and registered by the `postprocess` event, which was scheduled but had no handler. The figures were already registered by `plotFun()`, but the `plot` event called it without keeping the `sim` it returned, which discarded them. Tests cover the registration, the CSV bookkeeping and the `plot` event; the last fails if the result is dropped again.
* **Replaced the never-run test stub with metadata tests, and added testthat CI.** `tests/testthat/test-template.R` was the SpaDES boilerplate, unedited: paths from another machine, calls to `Event1`/`Event2` functions this module does not define, and assertions against placeholder strings. It had never been run and would have failed instantly, while making the module look tested. In its place, characterization tests over the module's public contract -- the input and output object names and classes, and the parameter names -- which is what a project binds to and what nothing checked until now. The expectations are GENERATED from the module's live metadata rather than transcribed, and were verified to fail when the contract changes.
* The contract test records that this module declares **no output objects**, which is by design rather than omission: everything it produces -- the per-refCode parquet aggregates, the envelope CSVs and the figures -- is written to `outputPath(sim)`, and nothing is handed back through the simList. A `createsOutput()` appearing later would be a real change in how it is consumed.

## Stable output keys via curated short names (`2.0.0.9020`)

- `refCode` (which names the `_aggregates/<refCode>/` parquet directory and the figure/CSV
  filenames) is now `LandWebUtils::refCodeFor(kind, layer)` = `<kind>_<slug(layer name)>`, instead
  of `paste0("<kind>_", abbreviate(layer, minlength = 8))`. `abbreviate()` was called one name at a
  time, and R only guarantees unique abbreviations *within a single call's vector* -- with the
  restored tenure x sub-region crossings, 68 of the 399 crossed FMA x ANSR names collided (e.g.
  `Canadian Forest Products Ltd. Alpine` and `Crowsnest Forest Products Ltd. Alpine` both ->
  `CFPLtd.A`), silently overwriting each other's aggregates and figures.
- Reporting-polygon layer names are now curated short names (`LandWebUtils::reportingPolygonLayers()$NAME_SHORT`,
  e.g. `ANSR`, or `SprayLake ANSR` for a crossed unit), so the on-disk key equals the label shown in
  figure titles and facet strips. **Existing `_aggregates/` and `figures/` directories carry the old
  keys and will be rebuilt.**
- Requires `LandWebUtils >= 1.0.3.9016`.

## Post-processing outputs reorganised; v2-form LandWeb-summary plots (`2.0.0.9009`)

- Multi-mode summary outputs now live under `outputs/<studyArea>/postprocess/` (a sibling of the
  mainSim rep dirs, which are still read from `outputPath(sim)`), with parallel `figures/` and `csv/`
  trees that mirror the same `<kind>/<layer>/` sub-structure (kind = lm / pm / boxplots / histograms;
  layer = the reporting-polygon-layer name). This replaces the flat `<refCode>_*.{csv,png}` files
  under `mainSim/`, so a figure and its data sit side by side and the dir path carries the context.
- The `lw` LandWeb summaries now render the v2 forms via `nrvtools::plot_leading_boxplot()` /
  `plot_largepatch_histogram()`: `figures/boxplots/<layer>/<subregion> <species>.png` (age-class
  boxplots + red CC dot) and `figures/histograms/<layer>/<size>/<subregion> <species>.png` (one file
  per species, four age-class panels, red CC line), replacing the busy distribution facets.
- `lm`/`pm` envelope figures relocate to `figures/<kind>/<layer>/{ribbon,boxplot}.png`.

## Save maps at `timeSeriesTimes` for the animation (`2.0.0.9008`)

- `mode = "single"` now generates and saves the `standAgeMap` / `vegTypeMap` at each
  `timeSeriesTimes` year (in addition to the summary times), so the `"am"` animation actually has
  its frames when read back in `mode = "multi"`. Previously the maps were only saved at
  `summaryPeriod` / `summaryInterval` times, so `timeSeriesTimes` outside that window produced no
  frames. The heavier `cohortData` / `pixelGroupMap` are still saved only at the summary analysis
  times, not at the (many) `timeSeriesTimes` years.

## Stand-age time-series animation (`2.0.0.9007`)

- New `"am"` `postprocessEvents` option (mode `"multi"`) ports the v2 `LandWeb_summary`
  animation: replicate 1's saved `standAgeMap` layers for the `timeSeriesTimes` years are
  reclassified into age classes (`ageClassCutOffs` / `ageClasses`, RdYlGn young &rarr; old,
  ggplot2 + `tidyterra`) and written to `standAge_animation.gif`. Encoded with **`gifski`**
  (pure-Rust) instead of `animation`/ImageMagick, so it needs no `policy.xml` cache
  configuration (the v2 path failed with "cache resources exhausted"; see LandWeb #153).
  Adds `gifski`, `RColorBrewer`, and `tidyterra` to `reqdPkgs`.

## `patchDirections` knob for lw large-patch connectivity (`2.0.0.9006`)

- New `patchDirections` parameter (integer, default `4L`) controls the patch
  connectivity of the `lw` large-patch analysis, passed through to
  `nrvtools::largePatchCounts()` / `landscapemetrics::get_patches()`: `4` = rook
  (4-connected) or `8` = queen (8-connected). **Departure from v2:** the v2
  `LandWeb_summary` large-patch analysis was fixed at 4-connectivity (GDAL
  `polygonize`'s default); the default `4L` reproduces v2, and `8L` (queen) is a
  deliberate v3-only option.

## LandWeb summaries: Leading + LargePatches (`2.0.0.9005`)

- `postprocess_lw` now runs the ported v2 `LandWeb_summary` analyses via `nrvtools` (>= 0.2.0):
  leading-vegetation-by-age-class and large-patch counts, keyed `lw_<poly>` (with a `_CC`
  current-conditions twin). The age basis is **time-since-fire** (`rstTimeSinceFire`, from
  burnSummaries), matching v2. The NRV distribution pools across replicates **and** summary years
  (summarized with `time` excluded from the id columns, via a new `.buildRepDataset(id_cols=)`
  passthrough). `plotFun` renders these as distribution histograms with a current-condition
  reference line (`plot_nrv_distribution()`) rather than time envelopes.

## Drop reporting-polygon features with no grouping label (`2.0.0.9003`)

* `landscapeMetrics` (lm) + `patchMetrics` (pm) now drop features whose grouping column (`Name`) is
  `NA` before summarising, and skip a layer entirely if none remain within the study area. An NA
  `polyName` made the metric producers select an empty subpoly (`summaryPolys[[col]] == NA`), so
  `crop()` returned NULL and the metric functions errored on `values(NULL)`.

## Landscape + patch metrics keyed on the reporting-polygon layer (`2.0.0.9002`)

* `landscapeMetrics` (lm) and `patchMetrics` (pm) now group by the `"Name"` column that
  `LandWebUtils::buildReportingPolygons()` sets and key each output `refCode` on the reporting-polygon
  **layer name** (`abbreviate(p)`), matching the `bc` (seral) event. Previously they hard-coded
  `polyCol = "NAME"` and keyed on `rptPoly[["ID"]]`, but the built reporting polygons carry neither a
  `NAME` nor an `ID` column, so a run errored inside `nrv_metrics_landscape`
  (`'names' attribute [1] must be the same length as the vector [0]`).

## Current-conditions reference = the saved year-0 state (`2.0.0.9001`)

* The `mode = "multi"` current-conditions (CC) snapshot now reads the simulation's saved year-0
  state (`rep01/vegTypeMap_year0000.tif` + `standAgeMap_year0000.tif`, the deterministic initial
  condition) instead of regenerating a VTM from `speciesLayers` in `patchMetrics` and reading a
  `reportingPolygons[["CC SAM"]]` stand-age. This removes (a) the write-before-read ordering bug
  where `landscapeMetrics` read the CC VTM before `patchMetrics` wrote it (a standalone run errored
  with `file does not exist: .../vegTypeMap_year0.tif`), and (b) the un-wired `"CC SAM"` dependency.
  The CC file paths are resolved once in `InitMulti`.

# NRV_summary 2.0.0

This is a breaking release that adopts the Arrow-native, memory-bounded NRV
summary path from `FOR-CAST/nrvtools (>= 0.1.0)` (which removed the former
in-memory `calculateLandscapeMetrics()` / `summarizePatchMetrics()` /
`summarizePatchMetricsSeral()`).

* landscape, patch, and seral patch metrics are now computed *raw* per replicate
  (`nrvtools::nrv_metrics_landscape()`, `calculatePatchMetrics()`,
  `calculatePatchMetricsSeral()`), written to per-replicate parquet partitions
  under `<outputPath>/_aggregates/<refCode>/replicate=<rep>/`, and reduced across
  replicates by `nrvtools::summarize_nrv()`, so the per-replicate rows are never
  all held in memory at once (the point of the change: large study areas x many
  replicates no longer blow up RAM).
* BREAKING output change: the former per-`funList` summary CSVs
  (`<refCode>_<funList>.csv`) and raw CSVs (`<refCode>_<funList>_raw.csv`) are
  replaced by a per-`refCode` envelope CSV (`<refCode>.csv`), one CSV per
  landscapemetrics metric (`<refCode>_<metric>.csv`), and the parquet dataset.
  The envelope columns are now `n_reps`/`mean`/`sd`/`min`/`q25`/`median`/`q75`/
  `max`/`se`/`ci` (the former `N`/`mn`/`mm`/`q1`/`md`/`q3`/`mx` are gone).
  Downstream readers must be updated.
* NRV plots use `nrvtools::plot_nrv_envelope()`, producing both a min-max ribbon
  (`<refCode>_ribbon.png`) and a box-and-whisker that shows the median and
  quartiles (`<refCode>_boxplot.png`) per reporting unit. The former paginated
  per-metric box/violin/over-time plots and the current-conditions overlay are
  pending re-addition (the `<refCode>_CC` envelope is still computed and saved).
* `SeralTable.csv` is unchanged in shape; its pooled total area is now
  `sum(n_reps * mean)` (identical to the former `sum(N * mn)`).
* removed the development `browser()` breakpoints from the postprocessing
  functions; requires `FOR-CAST/nrvtools (>= 0.1.0)`.

# NRV_summary 1.1.3

* current development version; consolidated NRV post-processing module that
  supersedes and absorbs the roles of `LandWeb_summary`, `timeSinceFire`, and
  `LandWeb_output`.
* run modes finalized: `mode = "single"` runs as part of a simulation (custom
  in-sim saving); `mode = "multi"` runs as post-processing across multiple runs.
* made `timeSeriesTimes` and `summaryPeriod` `start(sim)`-relative, with
  validation that both fall within the simulation's `start(sim)`..`end(sim)`
  range (reverted an interim change that had dropped the `start(sim)` offset).
* moved input-object existence checks out of `.inputObjects()` and into `init`.
* wired `standAgeMap` and `vegTypeMap` through `mod` and added them to the
  module's output objects; added missing module object dependencies.
* added `loadOrder` metadata; merged the upstream `startSim` branch;
  consistency tweaks.

# NRV_summary 1.1.2

* major refactor for combined single/multi use: added a `mode` parameter
  ("single" vs "multi") to select in-simulation vs post-processing behaviour,
  restructuring the module around the new params.
* replaced the `ml` map-list "god object" with a named `reportingPolygons`
  list; removed the `ml` output object and switched summaries to read polygon
  identity from the `reportingPolygons` "field" attribute.
* moved `summaryPeriod` from an expected input to a module parameter.
* added custom saving in single mode.

# NRV_summary 1.1.1

* added `sieveThresh` parameter and used `terra::sieve()` to merge singleton
  patches with neighbouring large patches; fixed `sieveThresh` usage and output
  filenames.
* bumped required `nrvtools` version for seral-stage map creation/fixes.
* removed the `sppEquivCol` parameter (cannot be user-defined).
* fixed grep regex used in file matching.

# NRV_summary 1.1.0

* minor release marking the stabilized seral-stage / `nrvtools`-based summary
  and mapping workflow.

# NRV_summary 1.0.0

* first stable release.
* seral-stage map creation and improved BC summaries against the latest
  `nrvtools`.
* fixed `refCodeCC` handling and `write.csv()` calls; removed redundant "CC"
  from CSV filenames; ungroup before `dplyr::summarize()` and pipe fixes;
  removed stray `browser()` calls.

# NRV_summary 0.0.6

* added `postprocessEvent` parameter to specify which post-processing events to
  run.
* added seral-stage patch figures; `standAgeMap` now taken from the `ml` object
  (`SAM`).
* stopped using package prefixes in `funList`.

# NRV_summary 0.0.1

* initial version: converted the BC NRV post-processing script into a SpaDES
  module (depends on `LandWebUtils`).
* landscape and patch metrics summaries by reporting polygon, including by LU,
  by BEC, and for arbitrary polygons supplied via `ml`; handled reporting
  subpolygons with empty data.frames and reported missing files.
* migrated to `terra`/`sf`; added violin plots and multi-page plots for patch
  metrics.
* parallel-processing improvements: core-count control via `future`, generalized
  `rasterListByPoly`, flammable map passed as a file for use in parallel; masked
  non-flammable pixels in TSF maps that otherwise skewed ages very old.
