#' Compute singular values of given matrix
#'
#' @description
#' `sing()` calculates the singular values of a given data matrix `D`. This is
#' done with a call to [svd()], and is included in `pcpr` to enable the quick
#' characterization of a data matrix's raw low-rank structure, to help decide
#' whether [rrmc()] or [root_pcp()] is the more appropriate PCP algorithm to
#' employ in conjunction with `D`.
#'
#' Experimentally, the [rrmc()] approach to PCP has best been able to handle
#' those datasets that are governed by complex underlying patterns characterized
#' by slowly decaying singular values, such as EH data. For observed data with a
#' well-defined low rank structure (rapidly decaying singular values),
#' [root_pcp()] may offer a better model estimate.
#'
#' @param D The input data matrix (cannot have `NA` values).
#'
#' @returns A numeric vector containing the singular values of `D`.
#' @examples
#' data <- sim_data()
#' sing(data$D)
#' # could plot the singular values for visual inspection with e.g.
#' # plot(sing(data$D), type = 'b')
#' @references "Singular value decomposition" [Wikipedia article](https://en.wikipedia.org/wiki/Singular_value_decomposition).
#' @export
sing <- function(D) {
  checkmate::assert_matrix(D, any.missing = FALSE)
  svd(D)$d
}
