pcp_rank_r = function(D,gamma,r) {

    # % solve the problem
    # %
    # %   min I_{rank(L) <=r} + gamma ||S||_1 + .5*||PO(Y-L-S)||_F^2
    # %
    # % using a proximal gradient method

    n = nrow(D)
    p = ncol(D)

    Om = !is.na(D)
    D[!Om] = 0

    L = matrix(0, nrow=n, ncol=p)
    S = matrix(0, nrow=n, ncol=p)

    obj = Inf

    t = .01
    iter = 0
    done = FALSE
    max_iter = 4000
    allObj = c()

    while (!done) {

        R = Om * ( D - L - S );
        S_new = prox_l1( S + t * R, t*gamma )
        L_new = proj_rank_r( L + t * R, r )

        delta = sqrt(norm(S-S_new,'F')^2 + norm(L-L_new,'F')^2)

        # norm 1 is the maximum absolute column sum of the matrix.
        obj_new = gamma * norm(S,"1") + .5 * norm( Om * ( D - L - S ), 'F' )^2


        if (obj_new > obj) {
            t = t * .95
        } else {
            t = t * 1.01
            S = S_new
            L = L_new
            delta_obj = obj - obj_new;
            obj = obj_new}

        iter = iter + 1

        if (iter %% 100 == 0) {print(paste('Iter:', iter, ', Obj:', obj, ', Step:', delta))}

        if (iter >= max_iter | t < 1e-10) {done = TRUE}

        allObj = rbind(allObj, obj)
    }
        return(list(L = L,S = S,iter = iter))
}




