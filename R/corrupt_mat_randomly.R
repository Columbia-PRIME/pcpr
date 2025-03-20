#' Corrupt given matrix with random missingness
#'
#' @description
#' `corrupt_mat_randomly()` corrupts a given data matrix `D` such that `perc`
#' percent of its entries are set to be missing (set to `NA`). Used by
#' [grid_search_cv()] in constructing test matrices for PCP models. Can be
#' used for experimentation with PCP models.
#'
#' Note: only _observed_ values can be corrupted as `NA`. This means if a matrix
#' `D` already has e.g. 20% of its values missing, then
#' `corrupt_mat_randomly(D, perc = 0.2)` would result in a matrix with 40% of
#' its values as missing.
#'
#' Should e.g. `perc = 0.6` be passed as input when `D` only has e.g. 10% of its
#' entries left as observed, then all remaining corruptable entries will be
#' set to `NA`.
#'
#' @param D The input data matrix.
#' @param perc A double `>= 0` specifying the percentage of entries in `D` to
#'   corrupt as missing (`NA`).
#' @param seed (Optional) An integer specifying the seed for the random
#'   selection of entries in `D` to corrupt as missing (`NA`). By default,
#'   `seed = 42`.
#'
#' @returns A list containing:
#' * `D_tilde`: The original matrix `D` with a random `perc` percent of its
#'   entries set to `NA`.
#' * `tilde_mask`: A binary matrix of `dim(D)` specifying the locations of
#'   corrupted entries (`1`) and uncorrupted entries (`0`).
#'
#' @seealso [grid_search_cv()], [sim_data()]
#' @examples
#' # Simple example corrupting 20% of a 5x5 matrix
#' D <- matrix(1:25, 5, 5)
#' corrupted_data <- corrupt_mat_randomly(D, perc = 0.2)
#' corrupted_data$D_tilde
#' sum(is.na(corrupted_data$D_tilde)) / prod(dim(corrupted_data$D_tilde))
#' # Now corrupting another 20% ontop of the original 20%
#' double_corrupted <- corrupt_mat_randomly(corrupted_data$D_tilde, perc = 0.2)
#' double_corrupted$D_tilde
#' sum(is.na(double_corrupted$D_tilde)) / prod(dim(double_corrupted$D_tilde))
#' # Corrupting the remaining entries by passing in a large value for perc
#' all_corrupted <- corrupt_mat_randomly(double_corrupted$D_tilde, perc = 1)
#' all_corrupted$D_tilde
#' @export
corrupt_mat_randomly <- function(D, perc, seed = 42) {
  # Calculate how many values need corrupting
  nvals_to_corrupt <- floor(perc * prod(dim(D)))
  D_vec <- as.vector(D)
  mask <- rep(0, length(D_vec))
  # Identify those values that can be corrupted
  pool <- which(!is.na(D_vec))
  if (length(pool) == 0) stop('There is nothing in the input matrix "D" that can be corrupted as missing.')
  # Corrupt values
  if (nvals_to_corrupt > length(pool)) {
    corrupted <- pool
  } else {
    set.seed(seed)
    corrupted <- sample(pool, nvals_to_corrupt, replace = FALSE)
  }
  mask[corrupted] <- 1
  D_vec[corrupted] <- NA
  # Reformat as matrix & return along w/binary mask specifying corruption locs
  n <- nrow(D)
  p <- ncol(D)
  ret_mat <- matrix(D_vec, n, p)
  ret_mask <- matrix(mask, n, p)
  list(D_tilde = ret_mat, tilde_mask = ret_mask)
}
