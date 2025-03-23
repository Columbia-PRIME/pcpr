#' Estimate sparsity of given matrix
#'
#' `sparsity()` estimates the percentage of entries in a given data matrix `D`
#' whose values are "practically zero". If the absolute value of an entry is
#' below a given threshold parameter `thresh`, then that value is determined
#' to be "practically zero", increasing the estimated sparsity of `D`. _Note
#' that `NA` values are imputed as 0 before the sparsity calculation is made._
#'
#' @param D The input data matrix.
#' @param thresh (Optional) A numeric threshold `>= 0` used to determine if an
#'   entry in `D` is "practically zero". If the absolute value of an entry is
#'   below `thresh`, then it is judged to be "practically zero". By default,
#'   `thresh = 1e-04`.
#'
#' @returns The sparsity of `D`, measured as the percentage of entries in `D`
#'   that are "practically zero".
#'
#' @seealso [matrix_rank()]
#' @examples
#' sparsity(matrix(rep(c(1, 0), 8), 4, 4))
#' sparsity(matrix(0:8, 3, 3))
#' sparsity(matrix(0, 3, 3))
#' @export
sparsity <- function(D, thresh = 1e-04) {
  checkmate::assert_matrix(D)
  checkmate::qassert(thresh, rules = "N1(0,)")
  if (sum(is.na(D)) > 0) D[is.na(D)] <- 0
  D[abs(D) < thresh] <- 0
  1 - (sum(D != 0) / prod(dim(D)))
}
