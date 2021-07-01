#' @export
RRMC = function(D, r, eta) {
	
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
	for (k in 1:r) {
		for (t in 0:T) {
			S = HT(Omega*(D - L), zeta)
			Mt = L + (n*p/sum(Omega)) * Omega * (D - S - L)
			L = proj_rank_r(Mt, k)
			temp = svd(Mt)$d
			zeta = eta * (temp[k+1] + 0.5^(t-2) * temp[k])
		}
		L_list[[k]] = L
		S_list[[k]] = S
	}

	list(L = L_list[[r]], S = S_list[[r]], L_list = L_list, S_list = S_list)
}