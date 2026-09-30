.onLoad <- function(libname, pkgname) {
  op <- options()
  op.latentcov.sbm <- list(
    latentcov.sbm.verbose = FALSE
  )

  # Only set options that are not already set by the user
  toset <- !(names(op.latentcov.sbm) %in% names(op))
  if (any(toset)) options(op.latentcov.sbm[toset])

  invisible()
}
