ilrInvcpp <- function(z, basis = default_Psi_function_cpp(ncol(z) + 1), log = FALSE) {
  return(ilrInv_cpp(z, basis, log = log))
}
