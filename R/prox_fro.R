# Prox Frobenius
prox_fro <- function(X,c) {

  n = norm(X,'F')

  if (n <= c) {Y = matrix(0, nrow = nrow(X), ncol = ncol(X))}
  else {Y = (1 - c/n) * X}

  return(Y)
}
