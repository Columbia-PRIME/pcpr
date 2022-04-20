#' @export
RRMC_block_missingness_guts = function(D, m1, m2, r, eta, numIter = 10, X11_L = NULL, verbose=FALSE) {
	
	n = nrow(D)
	p = ncol(D)

	Omega = !is.na(D)
	D[!Omega] = 0

	L_list = list()
	S_list = list()

	epsilon = 1e-05 * norm(D, "F")
	temp = svd(D)$d
	T = 10*log(40*r*n*temp[1]/epsilon)
	zeta = eta * temp[1]

	L = matrix(0, nrow = n, ncol = p)

	grid <- expand.grid(iter_t = 0:T, rank = 1:r)
	N <- nrow(grid)

	for (i in 1:N) {
		L_list[[i]] = L
		S = HT(Omega*(D - L), zeta)
		Mt = L + (n*p/sum(Omega)) * Omega * (D - S - L)
		L = proj_r_partial(Y=Mt, m1=m1, m2=m2, r=grid$rank[i], numIter=numIter, X11_L=X11_L, verbose=verbose)
		temp = svd(Mt)$d
		zeta = eta * (temp[grid$rank[i]+1] + 0.5^(grid$iter_t[i]-2) * temp[grid$rank[i]])
		cat(paste("Iter: ", i, " / ", N, "\n"))
	}
	L_list[[i]] = L

	list(L_list = L_list, grid = grid)

}