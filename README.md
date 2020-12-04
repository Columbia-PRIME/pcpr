# pcpr

PCP functions in R including adaptations for environmental data.

## includes:

1. [stable_pcp](R/stable_pcp.R)
    * This includes a non-negativity constraint on the `L` solution matrix. It does not include a LOD (limit of detection) penalty.
    * It takes 3 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter

2. [pcp_lod](R/pcp_lod.R)
    * This includes a non-negativity constraint on the `L` solution matrix and a separate penalty function for values <LOD.
    * Values <LOD should be pre-processed as `-1`.
    * It takes 3 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `LOD` may be a scalar, vector, or matrix
        
3. [root_pcp](R/root_pcp.R)
    * This changes the objective function by removing the squaring of the error term so that the `mu` parameter does not rely on the unknown `sigma` value. It does not include a non-negativity constraint on the `L` matrix or a LOD penalty.
    * It takes 3 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        
4. [root_pcp_nonnegL](R/root_pcp.R)
    * This includes a non-negativity constraint on the `L` matrix with the squareroot version of the objective function. It does not include a LOD penalty.
    * It takes 3 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
