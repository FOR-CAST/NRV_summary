Known issues: <https://github.com/FOR-CAST/NRV_summary/issues>

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
