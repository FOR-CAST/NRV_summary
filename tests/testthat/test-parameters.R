## vegLeadingProportion takes its default from LandR, as Biomass_core does: paramCheckOtherMods()
## stops a run whose modules disagree on it, and until 2.0.0.9028 this module's hard-coded 0.8
## disagreed with their 0.75.

paramDefault <- function(name) {
  p <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)$parameters
  p$default[[which(p$paramName == name)]]
}

test_that("vegLeadingProportion defaults to LandR::leadingSpeciesProp()", {
  expect_equal(paramDefault("vegLeadingProportion"), LandR::leadingSpeciesProp())
})

test_that("vegLeadingProportion follows the LandR.leadingSpeciesProp option", {
  withr::local_options(LandR.leadingSpeciesProp = 0.6)
  expect_equal(paramDefault("vegLeadingProportion"), 0.6)
})
