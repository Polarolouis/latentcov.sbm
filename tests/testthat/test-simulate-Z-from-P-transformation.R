# TDD: simulate_Z_from_P must honour the `transformation` argument.
# Bug: line 142 of R/simulators.R hardcodes `ilrInvcpp(P)` and ignores the
# provided transformation.

# a transformation that deterministically assigns every row to class 2
always_class_2 <- function(P) {
  matrix(c(0, 1, 0), nrow(P), 3, byrow = TRUE)
}

uniform_softmax <- function(P) {
  matrix(1 / 3, nrow(P), 3)
}

test_that("simulate_Z_from_P uses the provided transformation, not ilrInvcpp", {
  Z <- simulate_Z_from_P(P_fx, transformation = always_class_2)
  expect_identical(unique(as.integer(Z)), 2L)
})

test_that("simulate_Z_from_P honours the default ilrInvcpp transformation", {
  set.seed(601L)
  Z_default <- simulate_Z_from_P(P_fx)
  set.seed(601L)
  Z_explicit <- simulate_Z_from_P(P_fx, transformation = ilrInvcpp)
  expect_identical(as.integer(Z_default), as.integer(Z_explicit))
})
