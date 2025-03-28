
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pcpr <a href="https://columbia-prime.github.io/pcpr/"><img src="man/figures/logo.png" align="right" height="138" alt="pcpr website" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/pcpr)](https://cran.r-project.org/package=pcpr)
[![R-CMD-check](https://github.com/Columbia-PRIME/pcpr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Columbia-PRIME/pcpr/actions/workflows/R-CMD-check.yaml)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/pcpr)](https://cranlogs.r-pkg.org/badges/grand-total/pcpr)
<!-- badges: end -->

The R package `pcpr` implements Principal Component Pursuit (PCP), a
robust dimensionality reduction technique, for pattern recognition
tailored to environmental health data. The statistical methodology and
computational details are provided in Gibson et al. (2022).

## Installation

You can install the latest official CRAN release of `pcpr` with:

``` r
install.packages("pcpr")
```

The development version of `pcpr` can be installed from GitHub with:

``` r
# install.packages("pak")
pak::pak("Columbia-PRIME/pcpr")
```

`pcpr` can then be loaded and attached in your current R session as
usual with

``` r
library(pcpr)
```

## Getting help

Extensive documentation is available on our pkgdown
[website](https://columbia-prime.github.io/pcpr/reference/index.html)
and offline within R. You can see the `pcpr` reference manual in R with:

``` r
help("pcpr")
```

A number of vignettes are available from within R. They can be browsed
using:

``` r
browseVignettes("pcpr")
```

We recommend reading the vignettes in the following order:

1.  [Theory crash
    course](https://columbia-prime.github.io/pcpr/articles/theory-crash-course.html),
    or if directly in R: `vignette("theory-crash-course")`
2.  [Quickstart](https://columbia-prime.github.io/pcpr/articles/pcp-quickstart.html),
    or if directly in R: `vignette("pcp-quickstart")`
3.  [Air pollution source apportionment with
    PCP](https://columbia-prime.github.io/pcpr/articles/pcp-applied.html),
    or if directly in R: `vignette("pcp-applied")`

Have a bug to report or question to ask? [Open an issue on our
GitHub](https://github.com/Columbia-PRIME/pcpr/issues).

## Modeling overview

PCP algorithms model an observed exposure matrix $D$ as the sum of three
underlying ground-truth matrices:

![](man/figures/README-pcp-model-1.jpeg)

a low-rank matrix $L_0$ encoding consistent patterns of exposure, a
sparse matrix $S_0$ isolating unique or outlying exposure events (that
cannot be explained by the consistent exposure patterns), and dense
noise $Z_0$.

The models in `pcpr` seek to decompose an observed data matrix `D` into
estimated low-rank and sparse components `L` and `S` for use in
downstream environmental health analyses. The functions in `pcpr` are
outfitted with three environmental health (EH)-specific extensions
making `pcpr` particularly powerful for EH research:

1.  Missing value functionality
2.  Leveraging potential limit of detection (LOD) information
3.  Non-negativity constraint on the estimated `L` matrix

## PCP in environmental health studies

The methods in `pcpr` have already been applied in many environmental
health studies. Several are listed below:

- Tao et al. (2023) apply PCP to investigate the association between
  source-specific fine particulate matter and myocardial infarction
  hospitalizations in NYC.
- Wu et al. (2024) employ PCP for exposome profiling of environmental
  pollutants in seminal plasma, uncovering novel associations with semen
  parameters.
- Benavides et al. (2024) use PCP to develop a Community Severity Index
  in NYC, measuring the barrier effect of road infrastructure and
  traffic in cities.

## Acknowledgements

Please cite use of `pcpr` with:

Chillrud L, Benavides J, Gibson E, Zhang J, Yan J, Wright J, Goldsmith
J, Kioumourtzoglou M (2025). pcpr: Principal Component Pursuit for
Environmental Epidemiology. R package version 1.0.0,
<https://columbia-prime.github.io/pcpr/>,
<https://github.com/Columbia-PRIME/pcpr>.

    @Manual{,
      title = {pcpr: Principal Component Pursuit for Environmental Epidemiology},
      author = {Lawrence G. Chillrud and Jaime Benavides and Elizabeth A. Gibson and Junhui Zhang and Jingkai Yan and John N. Wright and Jeff Goldsmith and Marianthi-Anna Kioumourtzoglou},
      year = {2025},
      note = {R package version 1.0.0, https://github.com/Columbia-PRIME/pcpr},
      url = {https://columbia-prime.github.io/pcpr/},
    }

Please also cite Gibson et al. (2022).

This work was supported by NIEHS PRIME R01 ES028805.

Special thanks to Sophie Calhoun for designing `pcpr`’s logo!

## Usage

``` r
# In the below example, we simulate a simple mixtures model and run PCP,
# comparing it's performance to that of PCA. For an in depth example with
# simulated data, see vignette("pcp-quickstart"). For more realistic
# PCP usage, check out vignette("pcp-applied").

# Simulate an environmental mixture
data <- sim_data(
  n = 100, p = 10, r = 3,
  sparse_nonzero_idxs = seq(1, 1000, 101),
  sigma = 0.05
)
D <- data$D # Observed matrix
L_0 <- data$L # Ground truth low-rank matrix
S_0 <- data$S # Ground truth sparse matrix
Z_0 <- data$Z # Ground truth noise matrix

# Simulate a limit of detection for each chemical in mixture
lod_info <- sim_lod(D, q = 0.1)
D_lod <- lod_info$D_tilde
lod <- lod_info$lod

# Simulate missing observations
corrupted_data <- sim_na(D_lod, perc = 0.05)
D_tilde <- corrupted_data$D_tilde

# Finish simulating LOD by imputing values < LOD with LOD/sqrt(2)
lod_root2 <- matrix(
  lod / sqrt(2),
  nrow = nrow(D_tilde),
  ncol = ncol(D_tilde), byrow = TRUE
)
lod_idxs <- which(lod_info$tilde_mask == 1)
D_tilde[lod_idxs] <- lod_root2[lod_idxs]

# Run grid search to obtain optimal r, eta parameters
# (Not shown here to save space, see vignette("pcp-quickstart")
# for full example which obtains r = 3, eta = 0.224)
r_star <- 3
eta_star <- 0.224

# Run non-convex PCP to estimate L, S from D_tilde
pcp_model <- rrmc(D_tilde, r = r_star, eta = eta_star, LOD = lod)

# Clean up sparse matrix
pcp_model$S <- hard_threshold(pcp_model$S, thresh = 0.4)

# Benchmark with PCA's attempt at recovering L
D_imputed <- impute_matrix(D_tilde, apply(D_tilde, 2, mean, na.rm = TRUE))
L_pca <- proj_rank_r(D_imputed, r = r_star)

# Evaluate PCP ground truth
data.frame(
  "Obs_rel_err" = norm(L_0 - D_imputed, "F") / norm(L_0, "F"),
  "PCA_L_rel_err" = norm(L_0 - L_pca, "F") / norm(L_0, "F"),
  "PCP_L_rel_err" = norm(L_0 - pcp_model$L, "F") / norm(L_0, "F"),
  "PCP_S_rel_err" = norm(S_0 - pcp_model$S, "F") / norm(S_0, "F"),
  "PCP_L_rank" = matrix_rank(pcp_model$L),
  "PCP_S_sparsity" = sparsity(pcp_model$S)
)
#>   Obs_rel_err PCA_L_rel_err PCP_L_rel_err PCP_S_rel_err PCP_L_rank
#> 1   0.1440249    0.08096932    0.05847706      0.232115          3
#>   PCP_S_sparsity
#> 1          0.989
```

## References

Gibson, Elizabeth A., Junhui Zhang, Jingkai Yan, Lawrence Chillrud,
Jaime Benavides, Yanelli Nunez, Julie B. Herbstman, Jeff Goldsmith, John
Wright, and Marianthi-Anna Kioumourtzoglou. “Principal component pursuit
for pattern identification in environmental mixtures.” Environmental
Health Perspectives 130, no. 11 (2022): 117008.

Tao, Rachel H., Lawrence G. Chillrud, Yanelli Nunez, Sebastian T.
Rowland, Amelia K. Boehme, Jingkai Yan, Jeff Goldsmith, John Wright, and
Marianthi-Anna Kioumourtzoglou. “Applying principal component pursuit to
investigate the association between source-specific fine particulate
matter and myocardial infarction hospitalizations in New York City.”
Environmental Epidemiology 7 (2), (2023).

Wu, Haotian, Vrinda Kalia, Katherine E. Manz, Lawrence Chillrud,
Nathalie Hoffmann Dishon, Gabriela L. Jackson, Christian K. Dye, Raoul
Orvieto, Adva Aizer, Hagai Levine, Marianthi-Anna Kioumourtzoglou, Kurt
D. Pennell, Andrea A. Baccarelli, and Ronit Machtinger. “Exposome
Profiling of Environmental Pollutants in Seminal Plasma and Novel
Associations with Semen Parameters.” Environmental Science & Technology,
58 (31), (2024): 13594-13604.

Benavides, Jaime, Sabah Usmani, Vijay Kumar, and Marianthi-Anna
Kioumourtzoglou. “Development of a community severance index for urban
areas in the United States: A case study in New York City.” Environment
International, 185, (2024): 108526.
