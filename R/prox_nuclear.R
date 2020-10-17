#####
# Prox nuclear norm, L1 norm of the singular values
prox_nuclear <- function(Y,c) {

  USV <- svd(Y)
  U <- USV$u
  S <- USV$d
  V <- USV$v

  S_new <- sign(S) * pmax(abs(S) - c, 0)
  X <- U %*% diag(S_new) %*% t(V)
  nuclearX  <- sum(abs(S_new))

  return(list(X = X, nuclearX = nuclearX))
}
