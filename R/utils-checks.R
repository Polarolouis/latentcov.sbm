#' Check the emission distribution provided
#'
#' @param emission the string to check
#' @param arg the argument provided
#' @param call the calling context to go up the back trace
#' @keywords internal
check_emission_distribution <- function(
  emission,
  arg = rlang::caller_arg(emission),
  call = rlang::caller_env()
) {
  implemented_emissions <- c("poisson", "bernoulli")
  rlang::check_required(emission, arg = arg, call = call)
  rlang::arg_match(
    arg = emission,
    values = implemented_emissions,
    error_arg = arg,
    error_call = call
  )
}
