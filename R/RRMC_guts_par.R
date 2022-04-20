#' @importFrom magrittr %>%
#' @importFrom furrr future_map
#' @importFrom future plan multiprocess
RRMC_guts_par = function(D, r, eta) {
	
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
	
	future::plan(future::multiprocess)

	future_map(1:N, .init = L, .f = function(i) {
		S = HT(Omega*(D - L), zeta)
		Mt = L + (n*p/sum(Omega)) * Omega * (D - S - L)
		L = proj_rank_r(Mt, grid$rank[i])
		temp = svd(Mt)$d
		zeta = eta * (temp[grid$rank[i]+1] + 0.5^(grid$iter_t[i]-2) * temp[grid$rank[i]])
		cat(paste("Iter: ", i, " / ", N, "\n"))
	}) 

		L_list[[i]] = L
		S = HT(Omega*(D - L), zeta)
		Mt = L + (n*p/sum(Omega)) * Omega * (D - S - L)
		L = proj_rank_r(Mt, grid$rank[i])
		temp = svd(Mt)$d
		zeta = eta * (temp[grid$rank[i]+1] + 0.5^(grid$iter_t[i]-2) * temp[grid$rank[i]])
		cat(paste("Iter: ", i, " / ", N, "\n"))
	
	L_list[[i]] = L

	list(L_list = L_list, grid = grid)
}