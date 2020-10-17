# Prox L1 norm function, soft thresholding
prox_l1 <- function(Y, c) {
  X <- sign(Y) * pmax(abs(Y) - c, 0)
  return(X)
}
