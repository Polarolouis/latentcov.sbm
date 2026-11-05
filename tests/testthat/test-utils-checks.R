test_that("check_emission_distribution safely checks the emission", {
  # Detects non-existants distribution
  non_implem_distribs <- c("null", "chisq", "normal", "lognormal")
  for (distrib in non_implem_distribs) {
    expect_error(
      check_emission_distribution(distrib),
      regexp = glue::glue("not \"{distrib}\"")
    )
  }
  # Suggest corrections of typos
  typoed_distribs <- c("bernouilli", "bernouli", "poison", "poissson")
  names(typoed_distribs) <- c("bernoulli", "bernoulli", "poisson", "poisson")
  for (idx in seq_along(typoed_distribs)) {
    distrib <- typoed_distribs[idx]
    correct <- names(typoed_distribs)[idx]
    expect_error(
      check_emission_distribution(distrib),
      regexp = glue::glue("Did you mean \"{correct}\"?")
    )
  }
  # Works correctly for existing distribs
  implemented_emissions <- c("bernoulli", "poisson")
  for (distrib in implemented_emissions) {
    expect_no_error(check_emission_distribution(distrib))
    expect_equal(check_emission_distribution(distrib), distrib)
  }
})
