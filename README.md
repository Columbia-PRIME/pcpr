# pcpr

PCP functions in R including adaptations for environmental data.

## to run:

To use this package, clone the repo or download the .zip file. Locate the folder, unzip if applicable, and use the following code. Replace `path_to_folder` with your local path.

`install.packages("path_to_folder", repos = NULL, type="source")`   
`library(pcpr)`

## includes:

1. [stable_pcp](R/stable_pcp.R)
    * This includes a non-negativity constraint on the `L` solution matrix. It does not include a LOD (limit of detection) penalty.
    * It takes 4 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `verbose` parameter (optional)
    * For more info: `?stable_pcp`

2. [pcp_lod](R/pcp_lod.R)
    * This includes a non-negativity constraint on the `L` solution matrix and a separate penalty function for values <LOD.
    * Values <LOD should be pre-processed as `-1`.
    * It takes 5 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `LOD` may be a scalar, vector, or matrix
        * `verbose` parameter (optional)
    * For more info: `?pcp_lod`
        
3. [root_pcp](R/root_pcp.R)
    * This changes the objective function by removing the squaring of the error term so that the `mu` parameter does not rely on the unknown `sigma` value. It does not include a non-negativity constraint on the `L` matrix or a LOD penalty.
    * It takes 4 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `verbose` parameter (optional)
    * For more info: `?root_pcp`
        
4. [root_pcp_nonnegL](R/root_pcp_nonnegL.R)
    * This includes a non-negativity constraint on the `L` matrix with the squareroot version of the objective function. It does not include a LOD penalty.
    * It takes 4 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_nonnegL`
    
5. [root_pcp_na](R/root_pcp_na.R)
    * This allows for missing values with the squareroot version of the objective function. It does not include a LOD penalty or a non-negativity constraint on the `L` matrix.
    * Missing values should be pre-processed as `NA`.
    * It takes 4 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_na`
    
6. [root_pcp_na_nonnegL](R/root_pcp_na_nonnegL.R)
    * This includes a non-negativity constraint on the `L` matrix and allows for missing values with the squareroot version of the objective function. It does not include a LOD penalty.
    * Missing values should be pre-processed as `NA`.
    * It takes 4 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_na_nonnegL`
    
7. [root_pcp_na_nonnegL_lod](R/root_pcp_na_nonnegL_lod.R)
    * This includes a non-negativity constraint on the `L` matrix, allows for missing values, and includes a separate penalty function for values <LOD with the squareroot version of the objective function.
    * Missing values should be pre-processed as `NA`.
    * Values <LOD should be pre-processed as `-1`.
    * It takes 5 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `LOD` may be a scalar, vector, or matrix
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_na_nonnegL_lod`

8. [root_pcp_noncvx](R/root_pcp_noncvx.R)
    * This replaces the nuclear norm in the objective function with a projection to a lower rank. It does not include a LOD penalty or a non-negativity constraint on the `L` matrix.
    * It takes 5 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `r` the desired rank
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_noncvx`
    
9. [root_pcp_noncvx_nonneg](R/root_pcp_noncvx_nonneg.R)
    * This replaces the nuclear norm in the objective function with a projection to a lower rank. It includes a non-negativity constraint on the `L` matrix. It does not include a LOD penalty.
    * It takes 5 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `r` the desired rank
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_noncvx_nonneg`

10. [root_pcp_noncvx_w_na](R/root_pcp_noncvx_w_na.R)
    * This replaces the nuclear norm in the objective function with a projection to a lower rank. It does not include a LOD penalty or a non-negativity constraint on the `L` matrix. It does allow missing values.
    * Missing values should be pre-processed as `NA`.
    * It takes 5 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `r` the desired rank
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_noncvx_w_na`
    
11. [root_pcp_noncvx_nonneg_w_na](R/root_pcp_noncvx_nonneg_w_na.R)
    * This replaces the nuclear norm in the objective function with a projection to a lower rank. It includes a non-negativity constraint on the `L` matrix. It does not include a LOD penalty. It does allow missing values.
    * Missing values should be pre-processed as `NA`.
    * It takes 5 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `r` the desired rank
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_noncvx_nonneg_w_na`
    
12. [root_pcp_na_nonneg_noncvx_LOD](R/root_pcp_na_nonneg_noncvx_LOD.R)
    * This replaces the nuclear norm in the objective function with a projection to a lower rank. It includes a non-negativity constraint on the `L` matrix and a LOD penalty. It does allow missing values.
    * Missing values should be pre-processed as `NA`.
    * It takes 5 inputs:
        * `D` the original dataset
        * `lambda` parameter
        * `mu` parameter
        * `r` the desired rank
        * `LOD` may be a scalar, vector, or matrix
        * `verbose` parameter (optional)
    * For more info: `?root_pcp_na_nonneg_noncvx_LOD`
