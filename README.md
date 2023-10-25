# pcpr

PCP functions in R including adaptations for environmental data. This is the `dev` branch that [Lawrence](mailto:chili@u.northwestern.edu) is working on for PCP's software paper, and has been updated to run on R version 4.3.1 "Beagle Scouts". 

__Note: the `dev` branch may have bugs, incomplete documentation, etc. If you notice anything, please don't hestitate to submit an issue and tag @lawrence-chillrud, who will get to it as soon as possible. If you would like to make edits yourself, please do not edit dev directly. Instead, create an issue, start a new branch for that issue, and once that new branch is ready to be merged, request a code review. Once the review is approved, you are free to merge!__

## To install:

The simplest way to install `pcpr` is via the `devtools` package with:

```
install.packages("devtools") # if you don't already have devtools
devtools::install_github("https://github.com/Columbia-PRIME/pcpr/tree/dev") # you can also replace dev with any other branch you'd like to install
```

Once installed, you can load `pcpr` as you would any other package with `library(pcpr)`. 

## Main PCP functions:

Longer form documentation is coming soon, but some quick notes on the main functions in `pcpr`:

1. [stable_pcp](R/stable_pcp.R): Includes optional penalties for the limit of detection (LOD) penalty and a non-negativity constraint on the `L` solution matrix. Does not currently support missing values but will soon. It also needs documentation. Can be thought of as vanilla PCP with EHS enhancements (i.e., the optional penalties). See equation (15) in [Zhou et al. (2010)](https://ieeexplore.ieee.org/document/5513535) for the objective function of stable PCP, and look to [Gibson et al. (2022)](https://ehp.niehs.nih.gov/doi/pdf/10.1289/EHP10479) for details on the LOD and non-negativity penalties.

2. [root_pcp](R/root_pcp.R): Includes optional penalties for LOD and a non-negativity constraint on the `L` matrix. Accommodates missing values. Needs documentation. Can be thought of as the preferred method for data with a strong underlying low-rank structure characterized by rapidly decaying singular values (e.g., video or imaging data). See [Zhang, Yan, and Wright (2021)](https://proceedings.neurips.cc/paper/2021/hash/f65854da4622c1f1ad4ffeb361d7703c-Abstract.html) for the objective function of root PCP, and again [Gibson et al. (2022)](https://ehp.niehs.nih.gov/doi/pdf/10.1289/EHP10479) for the EHS-specific extensions. 
        
3. [RRMC](R/RRMC.R): LOD-penalty implemented, but the non-negativity constraint on the `L` matrix still needs to be written. Accommodates missingness. Needs documentation. Can be thought of as the preferred method for messy data that lacks a strong latent low-rank structure, characterized by singular values with a heavy tail (e.g., environmental mixtures). See [Cherapanamjeri, Gupta, and Jain (2017)](https://proceedings.mlr.press/v70/cherapanamjeri17a.html) for details and again [Gibson et al. (2022)](https://ehp.niehs.nih.gov/doi/pdf/10.1289/EHP10479) for the EHS-specific extensions.

I recommend looking at the [overview.Rmd](my-doc/overview.Rmd) file for a rough, unfinished, long-form documentation on how to use the `pcpr` package. There's also an associated [overview.html](my-doc/overview.html) that can be downloaded and should be able to be viewed in your browser for easier reading. Soon I'll incorporate a nicer version of this as a vignette, but for now hopefully the very rough html file will be better than nothing. You can also find some rough (and possibly outdated) documentation in the master branch for some of these functions. I'll flesh out the documentation in the coming weeks.

## Some utility functions:

1. [get_pcp_defaults](R/utils.R)

2. [grid_search_cv](R/grid_search_cv.R)

These should have some documentation but it is likely the help pages need updating - I'll get on this soon.