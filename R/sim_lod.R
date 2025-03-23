#' Simulate limit of detection data
#'
#' `sim_lod()` simulates putting the columns of a given matrix `D` under a limit
#' of detection (LOD) by calculating the given quantile `q` of each column and
#' corrupting all values < the quantile to `NA`, returning the newly corrupted
#' matrix, the binary corruption mask, and a vector of column LODs.
#'
#' @param D The input data matrix.
#' @param q A double in the range `[0, 1]` specifying the quantile to use
#'   in creating the column-wise LODs. Passed as the `probs` argument to the
#'   `quantile()` function.
#'
#' @returns A list containing:
#' * `D_tilde`: The original matrix `D` corrupted with < LOD `NA` values.
#' * `tilde_mask`: A binary matrix of `dim(D)` specifying the locations of
#'   corrupted entries (`1`) and uncorrupted entries (`0`).
#' * `lod`: A vector with `length(lod) == ncol(D)` providing the simulated
#'   LOD values corresponding to each column in the `D_tilde`.
#' @seealso [sim_na()], [impute_matrix()], [sim_data()]
#' @examples
#' D <- sim_data(5, 5, sigma = 0.8)$D
#' D
#' sim_lod(D, q = 0.2)
#' @export
#' @importFrom stats quantile
sim_lod <- function(D, q) {
  checkmate::assert_matrix(D)
  checkmate::qassert(q, "R1[0,1]")
  lod <- apply(D, MARGIN = 2, FUN = quantile, probs = q, na.rm = TRUE)
  lod_matrix <- matrix(lod, nrow = nrow(D), ncol = ncol(D), byrow = TRUE)
  lod_mask <- D < lod_matrix
  D[lod_mask] <- NA
  list(D_tilde = D, tilde_mask = 1 * lod_mask, lod = lod)
}
