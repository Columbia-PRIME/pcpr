# loss_lod function (only used in the loss function)
loss_lod <- function(X, D, LOD) {
  # % D is the original data
  # % X is the new thing (L + S)

  X_lod <- (X - D)   * (D >= 0) +
    (X - LOD)  * ((D < 0) & (X > LOD)) +
    X        * ((D < 0) & (X < 0))

  l <- sum(X_lod^2) / 2
  # % L2 norm

  # % Any D_ij < 0 AND X_ij < LOD AND > 0 are treated as equal
  # % Minimize discrepancy for valid data
  # % Want to shrink negative things

  l
}
