#' Cross-validated grid search for PCP models
#'
#' `grid_search_cv()` conducts a Monte Carlo style cross-validated grid search
#' of PCP parameters for a given data matrix `D`, PCP function `pcp_func`, and
#' grid of parameter settings to search through `grid`. The run time of the grid
#' search can be sped up using bespoke parallelization settings. See the below
#' sections for details.
#'
#' @section The Monte Carlo style cross-validation procedure:
#' Each hyperparameter setting is cross-validated by:
#' 1. Randomly corrupting `perc_test` percent of the entries in `D` as missing
#'   (i.e. `NA` values), yielding `D_tilde`. Done via [corrupt_mat_randomly()].
#' 2. Running the PCP function `pcp_func` on `D_tilde`, yielding estimates `L`
#'    and `S`.
#' 3. Recording the relative recovery error of `L` compared with the input data
#'    matrix `D` for _only those values that were imputed as missing during the
#'    corruption step_ (step 1 above). Mathematically, calculate:
#'    \eqn{||P_{\Omega^c}(D - L)||_F / ||P_{\Omega^c}(D)||_F}, where
#'    \eqn{P_{\Omega^c}} selects only those entries where
#'    `is.na(D_tilde) == TRUE`.
#' 4. Repeating steps 1-3 for a total of `runs` many times, where each "run" has
#'    a unique random seed from `1` to `runs` associated with it.
#' 5. Performance statistics can then be calculated for each "run", and then
#'    summarized across all runs for average model performance statistics.
#'
#' @section Best practices for `perc_test` and `runs`:
#' Experimentally, this grid search procedure retrieves the best performing
#' PCP parameter settings when `perc_test` is relatively low, e.g.
#' `perc_test = 0.05`, or 5%, and `runs` is relatively high, e.g. `runs = 100`.
#'
#' The larger `perc_test` is, the more the test set turns into a matrix
#' completion problem, rather than the desired matrix decomposition problem. To
#' better resemble the actual problem PCP will be faced with come inference
#' time, `perc_test` should therefore be kept relatively low.
#'
#' Choosing a reasonable value for `runs` is dependent on the need to keep
#' `perc_test` relatively low. Ideally, a large enough `runs` is used so that
#' many (if not all) of the entries in `D` are likely to eventually be tested.
#' Note that since test set entries are chosen randomly for all runs `1` through
#' `runs`, in the pathologically worst case scenario, the same exact test set
#' could be drawn each time. In the ebst case scenario, a different test set is
#' obtained each run, providing balanced coverage of `D`. Viewed another way,
#' the smaller `runs` is, the more the results are susceptible to overfitting to
#' the relatively few selected test sets.
#'
#' @section Interpretaion of results:
#' Once the grid search of has been conducted, the optimal hyperparameters can
#' be chosen by examining the output statistics `summary_stats`. Below are a
#' few suggestions for how to interpret the `summary_stats` table:
#' * Generally speaking, the first thing a user will want to inspect is the
#'   `rel_err` statistic, capturing the relative discrepancy between recovered
#'   test sets and their original, observed (yet possibly noisy) values. Lower
#'   `rel_err` means the PCP model was better able to recover the held-out test
#'   set. So, in general, **the best parameter settings are those with the
#'   lowest `rel_err`.** Having said this, it is important to remember that this
#'   statistic should be taken with a grain of salt: Because no ground truth `L`
#'   matrix exists, the `rel_err` measurement is forced to rely on the
#'   comparison between the _noisy observed data_ matrix `D` and the _estimated
#'   low-rank model_ `L`. So the `rel_err` metric is an "apples to oranges"
#'   relative error. For data that is a priori expected to be subject to a
#'   high degree of noise, it may actually be better to _discard_
#'   parameter settings with _suspiciously low_ `rel_err`s (in which
#'   case the solution may be hallucinating an inaccurate low-rank structure
#'   from the observed noise).
#' * For grid searches using [root_pcp()] as the PCP model, parameters that
#'   fail to converge can be discarded. Generally, fewer [root_pcp()] iterations
#'   (`num_iter`) taken to reach convergence portend a more reliable / stable
#'   solution. In rare cases, the user may need to increase [root_pcp()]'s
#'   `max_iter` argument to reach convergence. [RRMC()] does not report
#'   convergence metadata, as its optimization scheme runs for a fixed
#'   number of iterations.
#' * Parameter settings with unreasonable sparsity or rank measurements
#'   can also be discarded. Here, "unreasonable" means these reported metrics
#'   flagrantly contradict prior assumptions, knowledge, or work. For instance,
#'   most air pollution datasets contain a number of extreme exposure events, so
#'   PCP solutions returning sparse `S` models with 100% sparsity have obviously
#'   been regularized too heavily. Solutions with lower sparsities should be
#'   preferred. Note that reported sparsity and rank measurements are *estimates
#'   heavily dependent on the `thresh` set by the [sparsity()] & [matrix_rank()]
#'   functions*. E.g. it could be that the actual average matrix rank is much
#'   higher or lower when a threshold that better takes into account the
#'   relative scale of the singular values is used. Likewise for the sparsity
#'   estimations. Also, recall that the given value for `perc_test` artifically
#'   sets a sparsity floor, since those missing entries in the test set cannot
#'   be recovered in the `S` matrix. E.g. if `perc_test = 0.05`, then no
#'   parameter setting will have an estimated sparsity lower than 5%.
#'
#' @inheritParams rrmc
#' @param pcp_func The PCP function to use when grid searching. Must be either
#'   `rrmc` or `root_pcp` (passed without the soft brackets).
#' @param grid A `data.frame` of dimension `n` by `p` containing the `n`-many
#'   settings of `p`-many parameters to try. **The columns of `grid` should be
#'   named after the parameters in the function header of `pcp_func`.** For
#'   example, if `pcp_func = root_pcp`, then `names(grid)` must be set to
#'   `c("lambda", "mu")` if the grid is searching through different settings of
#'   `lambda` and `mu`. Likewise for `root_pcp = RRMC` and `eta` and `r`.
#' @param ... Any parameters required by `pcp_func` that could not be specified
#'   in `grid`. Importantly, these parameters are therefore kept constant (not
#'   involved in the grid search). The best example is the `LOD` parameter for
#'   those cases where the user has `LOD` information for PCP to leverage.
#' @param scale_func (Optional) The function used to scale the input `D` by
#'   column. By default, `scale_func = NULL`, and no scaling is to be done.
#' @param parallel_strategy (Optional) The parallelization strategy used when
#'   conducting the gridsearch (to be passed on to the [future::plan()]
#'   function). Must be one of: `"sequential"`, `"multisession"`, `"multicore"`
#'   or `"cluster"`. By default, `parallel_strategy = "multisession"`, which
#'   parallelizes the grid search via sockets in separate R _sessions_. If
#'   `parallel_strategy = "sequential"` then the search will be conducted in
#'   serial and the `cores` argument is ignored. The option
#'   `parallel_strategy = "multicore"` is not supported on Windows
#'   machines, nor in .Rmd files (must be run in a .R script) but parallelizes
#'   the search much faster than `"multisession"` since it runs separate
#'   _forked_ R processes. The option `parallel_strategy = "cluster"`
#'   parallelizes using separate R sessions running typically on one or more
#'   machines. Support for other parallel strategies will be added in a future
#'   release of `pcpr`. It is recommended to use
#'   `parallel_strategy = "multicore"` or `"multisession"` when possible.
#' @param cores (Optional) An integer specifying the number of cores to use when
#'   parallelizing the grid search. By default,
#'   `cores = parallel::detectCores(logical = F)`, which computes the number of
#'   physical CPUs available on the machine (see [parallel::detectCores()]).
#'   Ignored when `parallel_strategy = "sequential"`, must be `> 1` otherwise.
#' @param perc_test (Optional) The fraction of entries of `D` that will be
#'   randomly corrupted as `NA` missing values (the test set). Can be anthing in
#'   the range `[0, 1)`. By default, `perc_test = 0.05`. See **Best practices**
#'   section for more details.
#' @param runs (Optional) The number of times to test a given parameter setting.
#'   By default, `runs = 100`. See **Best practices** section for more details.
#' @param conserve_memory (Optional) A logical indicating if you only care about
#'   the statistics of the gridsearch and would therefore like to
#'   conserve memory when running the gridsearch. If set to `TRUE`, then only
#'   statistics on the parameters tested will be returned. By default,
#'   `conserve_memory = FALSE`, in which case additional objects saving the
#'   outputs of all runs of `pcp_func` will also be returned.
#' @param verbose (Optional) A logical indicating if you would like verbose
#'   output displayed or not (e.g. progress bars). By default, `verbose = TRUE`.
#' @param save_as (Optional) A character containing the root of the file path
#'   used to save the output to. Importantly, this should not end in any file
#'   extension, since this character will be used to save both the resulting
#'   `[save_as].rds` and `[save_as]_README.txt` files. By default,
#'   `save_as = NULL`, in which case the gridsearch is not saved to any file.
#'
#' @returns A list containing:
#' * `all_stats`: A `data.frame` containing the statistics of every run
#'   comprising the grid search. These statistics include the parameter
#'   settings for the run, along with the `run` number (used as the seed
#'   in the corruption step outlined in step 1 of the **Procedure** section),
#'   the relative error for the run `rel_err`, the rank of the recovered L
#'   matrix `L_rank`, the sparsity of the recovered S matrix `S_sparsity`,
#'   the number of `iterations` PCP took to reach convergence (for [root_pcp()]
#'   only), and the error status `run_error` of the PCP run (`NA` if no error,
#'   otherwise a character).
#' * `summary_stats`: A `data.frame` containing a summary of the information in
#'   `all_stats`. Summary made by column-wise averaging the results in
#'   `all_stats`.
#' * `L_mats`: A list containing all the `L` matrices returned from PCP
#'   throughout the gridsearch. Therefore, `length(L_mats) == nrow(all_stats)`.
#'   Row `i` in `all_stats` corresponds to `L_mats[[i]]`. Only returned when
#'   `conserve_memory = FALSE`.
#' * `S_mats`: A list containing all the S matrices returned from PCP throughout
#'   the gridsearch. Therefore, `length(S_mats) == nrow(all_stats)`. Row `i` in
#'   `all_stats` corresponds to `S_mats[[i]]`. Only returned when
#'   `conserve_memory = FALSE`.
#' * `test_mats`: A list of `length(runs)` containing all the corrupted test
#'   mats (and their masks) used throughout the gridsearch. Note:
#'   `all_stats$run[i]` corresponds to `test_mats[[i]]`. Only returned when
#'   `conserve_memory = FALSE`.
#' * `original_mat`: The original data matrix `D` _after it was column scaled by
#'   `scale_func`._ Only returned when `conserve_memory = FALSE`.
#' * `constant_params`: A copy of the constant parameters that were originally
#'   passed to the gridsearch (for record keeping).
#' @seealso [corrupt_mat_randomly()], [sparsity()], [matrix_rank()],
#'   [get_pcp_defaults()]
#' @examples
#' ####-------Simple simulated PCP problem-------####
#' # First we will simulate a simple dataset with the sim_data() function.
#' # The dataset will be a 100x10 matrix comprised of:
#' # 1. A rank-3 component as the ground truth L matrix;
#' # 2. A ground truth sparse component S w/outliers in 1st & last entries; and
#' # 3. A dense Gaussian noise component
#' data <- sim_data()
#' # Normally we would conduct grid search to tune eta. But, to keep the example
#' # short, we will just use best parameters from the below grid search example:
#' \dontrun{
#' eta_0 <- get_pcp_defaults(data$D)$eta
#' eta_grid <- data.frame("eta" = sort(c(0.1 * eta_0, eta_0 * seq(1, 10, 2))), "r" = 7)
#' gs <- grid_search_cv(data$D, rrmc, eta_grid)
#' dplyr::arrange(gs$summary_stats, rel_err)
#' }
#' # The gs found the best rank to be 3, and the best eta to be 0.3 or 0.4, so
#' # we will split the difference and use an eta of 0.35
#' pcp_model <- rrmc(data$D, r = 3, eta = 0.35)
#' data.frame(
#'   "Observed_relative_error" = norm(data$L - data$D, "F") / norm(data$L, "F"),
#'   "PCA_error" = norm(data$L - proj_rank_r(data$D, r = 3), "F") / norm(data$L, "F"),
#'   "PCP_L_error" = norm(data$L - pcp_model$L, "F") / norm(data$L, "F"),
#'   "PCP_S_error" = norm(data$S - pcp_model$S, "F") / norm(data$S, "F")
#' )
#' # Results:
#' # The grid search correctly found the rank (3) of the ground truth L matrix!
#' # PCP outperformed PCA in it's recovery of the L matrix (even though we let
#' # PCA "cheat" by telling PCA it was looking for a rank 3 solution)!
#' # PCP successfully isolated the outlying event in S!
#' @export
#' @importFrom magrittr %>%
grid_search_cv <- function(
  D,
  pcp_func,
  grid,
  ...,
  scale_func = NULL,
  parallel_strategy = "multisession",
  cores = parallel::detectCores(logical = FALSE),
  perc_test = 0.05,
  runs = 100,
  conserve_memory = FALSE,
  verbose = TRUE,
  save_as = NULL
) {
  # 0. Error handling:
  constant_params <- list(...)
  repeated_vars <- intersect(names(constant_params), colnames(grid))
  if (length(repeated_vars) > 0) stop(paste0('Arguments passed to "..." and "grid" are in conflict with one another. The following variables appear in both arguments and are therefore ambiguous: ', paste(repeated_vars, collapse = ", ")))
  if (!is.matrix(D)) stop('Invalid value passed for argument "D". Argument "D" must be a matrix. See documentation with "?grid_search_cv".')
  if (!(parallel_strategy %in% c("sequential", "multisession", "multicore", "cluster"))) stop('Invalid value passed for argument "parallel_strategy". Must be one of: {"sequential", "multisession", "multicore", "cluster"}. See documentation with "?grid_search_cv".')
  if (parallel_strategy != "sequential" && cores == 1) stop('Arguments "parallel_strategy" and "cores" are in conflict with one another. If you want to use 1 core then "parallel_strategy" should be set to "sequential". If you want to run the gridsearch in parallel then cores should be > 1.')
  if (parallel_strategy == "multicore" && Sys.info()['sysname'] == "Windows") stop('Argument "parallel_strategy" cannot be set to "multicore" on a Windows machine. Please use "multisession", "cluster", or "sequential" instead. See documentation with "?grid_search_cv".')
  if (perc_test < 0 || perc_test >= 1) stop('Invalid value passed for argument "perc_test". Argument "perc_test" must be in the range [0, 1). See documentation with "?grid_search_cv".')
  if (runs < 1) stop('Invalid value passed for argument "runs". Argument "runs" must be >= 1. See documentation with "?grid_search_cv".')
  if (perc_test == 0 && runs > 1) stop('Arguments "perc_test" and "runs" are in conflict with one another. Argument "runs" cannot be greater than 1 if "perc_test" is 0. See documentation with "?grid_search_cv".')
  if (!is.logical(verbose)) stop('Invalid value passed for argument "verbose". Must be a logical (TRUE / FALSE). See documentation with "?grid_search_cv".')
  if (!is.null(save_as) && stringr::str_detect(save_as, "\\.")) stop('Invalid value passed for argument "save_as". There must NOT be any periods (".") in "save_as". The proper file extensions are automatically added.')
  # 1. Setting up the search:
  if (verbose) cat("\nInitializing gridsearch...")
  # 1a. File names for saving the search later:
  if (!is.null(save_as)) {
    README_file <- paste0(save_as, "_README.txt")
    RDS_file <- paste0(save_as, ".rds")
    file_version <- 1
    while (file.exists(README_file) || file.exists(RDS_file)) {
      README_file <- paste0(save_as, file_version, "_README.txt")
      RDS_file <- paste0(save_as, file_version, ".rds")
      file_version <- file_version + 1
    }
    if (verbose) cat(paste0("\nThe completed gridsearch will be saved to the following files:\n\t", RDS_file, "\n\t", README_file))
  } else if (verbose) cat("\nThe completed gridsearch will NOT be saved to any files, but simply returned.")
  # 1b. Parsing the grid that was passed:
  metrics <- c("rel_err", "L_rank", "S_sparsity", "iterations", "run_error", "run_error_perc")
  for (metric in metrics) {
    if (!metric %in% names(grid)) {
      grid[, metric] <- NA
    }
  }
  # 1c. Setting up the params to search through:
  param_names <- grid %>% dplyr::select(!tidyselect::all_of(metrics)) %>% colnames()
  points_to_eval <- which(is.na(grid$rel_err))
  params <- data.frame(grid[points_to_eval, param_names], run_num = rep(1:runs, each = length(points_to_eval)), row.names = NULL)
  colnames(params) <- c(param_names, "run_num")
  # 1d. Scaling the matrix by column:
  if (!is.null(scale_func)) D <- apply(D, 2, scale_func)
  # 1e. Creating the test sets:
  test_mats <- purrr::map(1:runs, ~ corrupt_mat_randomly(D, perc = perc_test, seed = .x))
  # 1f. Parallel programming setup:
  start_time <- Sys.time()
  if (cores == 1) {
    future::plan(parallel_strategy)
    if (verbose) cat(paste0("\nBeginning sequential gridsearch...\nStart time: ", start_time, "\n"))
  } else {
    future::plan(parallel_strategy, workers = cores)
    if (verbose) cat(paste0("\nBeginning parallel gridsearch using ", cores, " cores and a ", parallel_strategy, " strategy...\nStart time: ", start_time, "\n"))
  }
  # 1g. Progress bar setup:
  p <- progressr::progressor(steps = nrow(params))
  old_progress_bar <- verbose && parallel_strategy == "multicore"
  # 2. Conducting the search:
  if (conserve_memory) {
    pcp_evals <- furrr::future_map_dfr(1:nrow(params), .f = function(i) {
      p()
      tryCatch(expr = {
        param_setting <- as.data.frame(params[i, param_names])
        colnames(param_setting) <- param_names
        pcp_mat <- do.call(pcp_func, c(constant_params, as.list(param_setting), list(D = test_mats[[params$run_num[i]]]$D_tilde)))
        eval_params(settings = params[i,], test_mat = D, pcp_model = pcp_mat, test_mask = test_mats[[params$run_num[i]]]$tilde_mask)
      }, error = function(e) {
        cbind(params[i,], data.frame(rel_err = NA, L_rank = NA, S_sparsity = NA, iterations = NA, run_error = paste0(e)))
      })
    }, .progress = old_progress_bar)
    future::plan("sequential")
  } else {
    pcp_mats <- furrr::future_map(1:nrow(params), .f = function(i) {
      p()
      param_setting <- as.data.frame(params[i, param_names])
      colnames(param_setting) <- param_names
      tryCatch(expr = {
        do.call(pcp_func, c(constant_params, as.list(param_setting), list(D = test_mats[[params$run_num[i]]]$D_tilde)))
      }, error = function(e) {
        list(message = paste0(e), settings = as.list(param_setting), test_mat = test_mats[[params$run_num[i]]]$D_tilde)
      })
    }, .progress = old_progress_bar)
    future::plan("sequential")
    pcp_evals <- purrr::map_dfr(1:nrow(params), .f = function(i) {
      eval_params(settings = params[i,], test_mat = D, pcp_model = pcp_mats[[i]], test_mask = test_mats[[params$run_num[i]]]$tilde_mask)
    })
  }
  end_time <- Sys.time()
  if (verbose) cat(paste0("\nGridsearch completed at time: ", end_time))
  if (perc_test == 0) pcp_evals$rel_err <- NA
  # 3. Summarizing the results from the search:
  if (runs > 1) {
    evals_summary <- pcp_evals %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(param_names))) %>%
      dplyr::summarise(rel_err = mean(rel_err, na.rm = T), L_rank = mean(L_rank, na.rm = T), S_sparsity = mean(S_sparsity, na.rm = T), iterations = mean(iterations, na.rm = T), run_error_perc = paste0(round(100 * sum(!is.na(run_error)) / dplyr::n(), 2), "%"), .groups = "drop")
  } else {
    evals_summary <- pcp_evals %>% dplyr::select(!run_num)
  }
  if (verbose) cat("\nMetrics calculations complete.")
  # 4. package results:
  if (conserve_memory) {
    results <- list(all_stats = pcp_evals, summary_stats = evals_summary, constant_params = constant_params)
  } else {
    if ("L_list" %in% names(pcp_mats[[1]])) {
      Ls <- foreach::`%do%`(
        foreach::foreach(k = 1:length(pcp_mats), .combine = c),
        pcp_mats[[k]]$L_list
      )
      Ss <- foreach::`%do%`(
        foreach::foreach(k = 1:length(pcp_mats), .combine = c),
        pcp_mats[[k]]$S_list
      )
    } else {
      Ls <- foreach::`%do%`(
        foreach::foreach(k = 1:length(pcp_mats), .combine = c),
        pcp_mats[[k]]$L
      )
      Ss <- foreach::`%do%`(
        foreach::foreach(k = 1:length(pcp_mats), .combine = c),
        pcp_mats[[k]]$S
      )
    }
    results <- list(all_stats = pcp_evals, summary_stats = evals_summary, L_mats = Ls, S_mats = Ss, test_mats = test_mats, original_mat = D, constant_params = constant_params)
  }
  # 5. save results & write README file:
  if (!is.null(save_as)) {
    if (verbose) cat(paste0("\nNow saving the completed gridsearch to the following files:\n\t", RDS_file, "\n\t", README_file))
    saveRDS(results, file = RDS_file)
    params_summary <- c()
    for (p in param_names) params_summary <- c(params_summary, paste0("\t\t", p, ": {", paste(sort(unique(params[[p]])), collapse = ", "), "}"))
    README <- file(README_file)
    writeLines(c(
      paste(rep('-', 80), collapse =""),
      paste0("\nFile: ", README_file),
      paste0("Date & Time: ", Sys.time()),
      paste0("Description: README file for the gridsearch saved to ", RDS_file, "\n"),
      paste(rep('-', 80), collapse =""),
      paste0("\tGridsearch start time: ", start_time),
      paste0("\tGridsearch end time: ", end_time),
      paste0("\tPCP function used: ", deparse(pcp_func)[1]),
      paste0("\tTotal number of settings searched through: ", nrow(params)),
      paste0("\tParameters searched through: "),
      params_summary,
      paste0("\nGridsearch settings:"),
      paste0("\tScale function used: ", deparse(scale_func)[1]),
      paste0("\tParallelization strategy used: ", parallel_strategy),
      paste0("\tCores used: ", cores),
      paste0("\tPercent of matrix corrupted for test set: ", round(perc_test*100, 3), "%"),
      paste0("\tRuns per parameter setting: ", runs),
      paste0("\tMemory conserved?: ", conserve_memory),
      paste0("\tConstant parameters passed to PCP: ", paste(names(constant_params), collapse = ", "), "\n"),
      paste(rep("-", 80), collapse = "")
    ), README)
    close(README)
    if (verbose) cat("\nSave completed. Exiting gridsearch.")
  }
  # 6. Return results:
  results
}

#' Corrupt given matrix with random missingness
#'
#' @description
#' `corrupt_mat_randomly()` corrupts a given data matrix `D` such that `perc`
#' percent of its entries are set to be missing (set to `NA`). Used by
#' [grid_search_cv()] in constructing test matrices for PCP models. Can be
#' used for experimentation with PCP models.
#'
#' Note: only _observed_ values can be corrupted as `NA`. This means if a matrix
#' `D` already has e.g. 20% of its values missing, then
#' `corrupt_mat_randomly(D, perc = 0.2)` would result in a matrix with 40% of
#' its values as missing.
#'
#' Should e.g. `perc = 0.6` be passed as input when `D` only has e.g. 10% of its
#' entries left as observed, then all remaining corruptable entries will be
#' set to `NA`.
#'
#' @param D The input data matrix.
#' @param perc A double `>= 0` specifying the percentage of entries in `D` to
#'   corrupt as missing (`NA`).
#' @param seed (Optional) An integer specifying the seed for the random
#'   selection of entries in `D` to corrupt as missing (`NA`). By default,
#'   `seed = 42`.
#'
#' @returns A list containing:
#' * `D_tilde`: The original matrix `D` with a random `perc` percent of its
#'   entries set to `NA`.
#' * `tilde_mask`: A binary matrix of `dim(D)` specifying the locations of
#'   corrupted entries (`1`) and uncorrupted entries (`0`).
#'
#' @seealso [grid_search_cv()], [sim_data()]
#' @examples
#' # Simple example corrupting 20% of a 5x5 matrix
#' D <- matrix(1:25, 5, 5)
#' corrupted_data <- corrupt_mat_randomly(D, perc = 0.2)
#' corrupted_data$D_tilde
#' sum(is.na(corrupted_data$D_tilde)) / prod(dim(corrupted_data$D_tilde))
#' # Now corrupting another 20% ontop of the original 20%
#' double_corrupted <- corrupt_mat_randomly(corrupted_data$D_tilde, perc = 0.2)
#' double_corrupted$D_tilde
#' sum(is.na(double_corrupted$D_tilde)) / prod(dim(double_corrupted$D_tilde))
#' # Corrupting the remaining entries by passing in a large value for perc
#' all_corrupted <- corrupt_mat_randomly(double_corrupted$D_tilde, perc = 1)
#' all_corrupted$D_tilde
#' @export
corrupt_mat_randomly <- function(D, perc, seed = 42) {
  # Calculate how many values need corrupting
  nvals_to_corrupt <- floor(perc * prod(dim(D)))
  D_vec <- as.vector(D)
  mask <- rep(0, length(D_vec))
  # Identify those values that can be corrupted
  pool <- which(!is.na(D_vec))
  if (length(pool) == 0) stop('There is nothing in the input matrix "D" that can be corrupted as missing.')
  # Corrupt values
  if (nvals_to_corrupt > length(pool)) {
    corrupted <- pool
  } else {
    set.seed(seed)
    corrupted <- sample(pool, nvals_to_corrupt, replace = FALSE)
  }
  mask[corrupted] <- 1
  D_vec[corrupted] <- NA
  # Reformat as matrix & return along w/binary mask specifying corruption locs
  n <- nrow(D)
  p <- ncol(D)
  ret_mat <- matrix(D_vec, n, p)
  ret_mask <- matrix(mask, n, p)
  list(D_tilde = ret_mat, tilde_mask = ret_mask)
}

#' Evaluate parameters
#'
#' Given documented parameter settings `settings`, a test matrix `test_mat`,
#' PCP solution `pcp_model`, and binary mask `test_mask`, returns statistics
#' about the `pcp_model`'s performance with `settings` on the `test_mat`.
#' **This is an internal function needed by [grid_search_cv()]. It is not
#' expected that users should require access to this function.**
#'
#' @param settings Documented parameter settings associated with the given
#'   `pcp_model`.
#' @param test_mat The given test matrix with `test_mask` _observed_ (i.e. not
#'   missing).
#' @param pcp_model The output of a PCP algorithm, either [rrmc()] or
#'   [root_pcp()].
#' @param test_mask The binary mask indicating which entries in `test_mat`
#'   comprise the test set.
#'
#' @returns A `data.frame` object with PCP's performance statistics.
#'
#' @seealso [grid_search_cv()]
#' @keywords internal
eval_params <- function(settings, test_mat, pcp_model, test_mask) {
  if ("message" == names(pcp_model)[1]) {
    stats <- cbind(settings, data.frame(rel_err = NA, L_rank = NA, S_sparsity = NA, iterations = NA, run_error = pcp_model$message))
  } else if ("L_list" == names(pcp_model)[3]) {
    if ("r" %in% names(settings)) settings <- settings[setdiff(names(settings), "r")]
    test_mat[is.na(test_mat)] <- 0
    stats <- purrr::imap_dfr(.x = seq(1, length(pcp_model$L_list)), ~ data.frame(
      r = .y, settings,
      rel_err = norm((test_mat - pcp_model$L_list[[.x]])*test_mask, "F") / norm(test_mat*test_mask, "F"),
      L_rank = matrix_rank(pcp_model$L_list[[.x]], 1e-04),
      S_sparsity = sparsity(pcp_model$S_list[[.x]]),
      iterations = NA, run_error = NA
    ))
  } else {
    test_mat[is.na(test_mat)] <- 0
    stats <- settings
    stats$rel_err <- norm((test_mat - pcp_model$L)*test_mask, "F") / norm(test_mat*test_mask, "F")
    stats$L_rank <- matrix_rank(pcp_model$L, 1e-04)
    stats$S_sparsity <- sparsity(pcp_model$S)
    stats$iterations <- pcp_model$num_iter
    stats$run_error <- NA
  }
  stats
}

#' Estimate sparsity of given matrix
#'
#' `sparsity()` estimates the percentage of entries in a given data matrix `D`
#' whose values are "practically zero". If the absolute value of an entry is
#' below a given threshold parameter `thresh`, then that value is determined
#' to be "practically zero", increasing the estimated sparsity of `D`.
#'
#' @param D The input data matrix.
#' @param thresh (Optional) A numeric threshold `>= 0` used to determine if an
#'   entry in `D` is "practically zero". If the absolute value of an entry is
#'   below `thresh`, then it is judged to be "practically zero". By default,
#'   `thresh = 1e-04`.
#'
#' @returns The sparsity of `D`, measured as the percentage of entries in `D`
#'   that are "practically zero".
#'
#' @seealso [rank_matrix()]
#' @examples
#' sparsity(matrix(rep(c(1, 0), 8), 4, 4))
#' sparsity(matrix(0:8, 3, 3))
#' sparsity(matrix(0, 3, 3))
#' @export
sparsity <- function(D, thresh = 1e-04) {
  D[abs(D) < thresh] <- 0
  100 - (100 * sum(D != 0) / prod(dim(D)))
}

#' Estimate rank of a given matrix
#'
#' @description
#' `matrix_rank()` estimates the rank of a given data matrix `D` by counting the
#' number of "practically nonzero" singular values of `D`.
#'
#' A singular value \eqn{s} is determined to be "practically nonzero" if
#' \eqn{s \geq s_max * thresh}, i.e. if it is greater than or equal to the
#' maximum singular value in `D` scaled by a given threshold `thresh`.
#'
#' @param D The input data matrix.
#' @param thresh (Optional) A double `eqn{> 0}`, specifying the relative
#'   threshold by which "practically zero" is determined, used to calculate the
#'   rank of `D`. By default, `thresh = NULL`, in which case the threshold is
#'   set to `max(dim(D)) * .Machine$double.eps`.
#'
#' @returns An integer representing the rank of `D`.
#'
#' @seealso [sparsity()]
#' @examples
#' data <- sim_data()
#' matrix_rank(data$D)
#' matrix_rank(data$L)
#' @export
matrix_rank <- function(D, thresh = NULL) {
  if (is.null(thresh)) thresh <- max(dim(D)) * .Machine$double.eps
  singular_values <- svd(D)$d
  sum(singular_values >= max(singular_values) * thresh)
}