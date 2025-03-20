#' Estimate rank of a given matrix
#'
#' @description
#' `matrix_rank()` estimates the rank of a given data matrix `D` by counting the
#' number of "practically nonzero" singular values of `D`.
#'
#' The rank of a matrix is the number of linearly independent columns or rows in
#' the matrix, governing the structure of the data. It can intuitively be
#' thought of as the number of inherent latent patterns in the data.
#'
#' A singular value \eqn{s} is determined to be "practically nonzero" if
#' \eqn{s \geq s_{max} \cdot thresh}, i.e. if it is greater than or equal to the
#' maximum singular value in `D` scaled by a given threshold `thresh`.
#'
#' @param D The input data matrix.
#' @param thresh (Optional) A double \eqn{> 0}, specifying the relative
#'   threshold by which "practically zero" is determined, used to calculate the
#'   rank of `D`. By default, `thresh = NULL`, in which case the threshold is
#'   set to `max(dim(D)) * .Machine$double.eps`.
#'
#' @returns An integer estimating the rank of `D`.
#'
#' @seealso [sparsity()]
#' @examples
#' data <- sim_data()
#' matrix_rank(data$D)
#' matrix_rank(data$L)
#' @export
matrix_rank <- function(D, thresh = NULL) {
  if (is.null(thresh)) thresh <- max(dim(D)) * .Machine$double.eps
  singular_values <- svd(D)$d
  sum(singular_values >= max(singular_values) * thresh)
}
