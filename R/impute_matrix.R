#' Impute missing values in given matrix
#'
#' `impute_matrix()` imputes the missing `NA` values in a given matrix using a
#' given `imputation_scheme`.
#'
#' @param D The input data matrix.
#' @param imputation_scheme The values to replace missing `NA` values in `D`
#'   with. Can be either:
#'   * A scalar numeric, indicating all `NA` values should be imputed with the
#'     same scalar numeric value;
#'   * A vector of length `ncol(D)`, signifying column-specific imputation,
#'     where each entry in the `imputation_scheme` vector corresponds to the
#'     imputation value for each column in `D`; or
#'   * A matrix of dimension `dim(D)`, indicating an observation-specific
#'     imputation scheme, where each entry in the `imputation_scheme` matrix
#'     corresponds to the imputation value for each entry in `D`.
#'
#' @returns The imputed matrix.
#' @seealso [sim_na()], [sim_lod()], [sim_data()]
#' @examples
#' #### ------------Imputation with a scalar------------####
#' # simulate a small 5x5 mixture
#' D <- sim_data(5, 5)$D
#' # corrupt the mixture with 40% missing observations
#' D_tilde <- sim_na(D, 0.4)$D_tilde
#' D_tilde
#' # impute missing values with 0
#' impute_matrix(D_tilde, 0)
#' # impute missing values with -1
#' impute_matrix(D_tilde, -1)
#'
#' #### ------------Imputation with a vector------------####
#' # impute missing values with the column-mean
#' impute_matrix(D_tilde, apply(D_tilde, 2, mean, na.rm = TRUE))
#' # impute missing values with the column-min
#' impute_matrix(D_tilde, apply(D_tilde, 2, min, na.rm = TRUE))
#'
#' #### ------------Imputation with a matrix------------####
#' # impute missing values with random Gaussian noise
#' noise <- matrix(rnorm(prod(dim(D_tilde))), nrow(D_tilde), ncol(D_tilde))
#' impute_matrix(D_tilde, noise)
#'
#' #### ------------Imputation with LOD/sqrt(2)------------####
#' D <- sim_data(5, 5)$D
#' lod_info <- sim_lod(D, q = 0.2)
#' D_tilde <- lod_info$D_tilde
#' D_tilde
#' lod <- lod_info$lod
#' impute_matrix(D_tilde, lod / sqrt(2))
#' @export
impute_matrix <- function(D, imputation_scheme) {
  checkmate::assert_matrix(D)
  checkmate::qassert(imputation_scheme, rules = "n+")
  if (checkmate::test_numeric(imputation_scheme, min.len = 2) && !checkmate::test_matrix(imputation_scheme)) {
    checkmate::assert_true(ncol(D) == length(imputation_scheme))
  }
  if (checkmate::test_matrix(imputation_scheme)) {
    checkmate::assert_true(all(dim(D) == dim(imputation_scheme)))
  }
  n <- nrow(D)
  p <- ncol(D)
  if (!checkmate::test_matrix(imputation_scheme)) {
    imputation_scheme <- matrix(imputation_scheme, nrow = n, ncol = p, byrow = TRUE)
  }
  D[is.na(D)] <- imputation_scheme[is.na(D)]
  D
}
