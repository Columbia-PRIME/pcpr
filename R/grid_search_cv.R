#' Cross-validated grid search for PCP models
#'
#' `grid_search_cv()` conducts a Monte Carlo style cross-validated grid search
#' of PCP parameters for a given data matrix `D`, PCP function `pcp_fn`, and
#' grid of parameter settings to search through `grid`. The run time of the grid
#' search can be sped up using bespoke parallelization settings. The call to
#' `grid_search_cv()` can be wrapped in a call to [progressr::with_progress()]
#' for progress bar updates. See the below sections for details.
#'
#' @section The Monte Carlo style cross-validation procedure:
#' Each hyperparameter setting is cross-validated by:
#' 1. Randomly corrupting `perc_test` percent of the entries in `D` as missing
#'   (i.e. `NA` values), yielding `D_tilde`. Done via [sim_na()].
#' 2. Running the PCP function `pcp_fn` on `D_tilde`, yielding estimates `L`
#'    and `S`.
#' 3. Recording the relative recovery error of `L` compared with the input data
#'    matrix `D` for _only those values that were imputed as missing during the
#'    corruption step_ (step 1 above). Mathematically, calculate:
#'    \eqn{||P_{\Omega^c}(D - L)||_F / ||P_{\Omega^c}(D)||_F}, where
#'    \eqn{P_{\Omega^c}} selects only those entries where
#'    `is.na(D_tilde) == TRUE`.
#' 4. Repeating steps 1-3 for a total of `num_runs`-many times, where each "run"
#'    has a unique random seed from `1` to `num_runs` associated with it.
#' 5. Performance statistics can then be calculated for each "run", and then
#'    summarized across all runs for average model performance statistics.
#'
#' @section Best practices for `perc_test` and `num_runs`:
#' Experimentally, this grid search procedure retrieves the best performing
#' PCP parameter settings when `perc_test` is relatively low, e.g.
#' `perc_test = 0.05`, or 5%, and `num_runs` is relatively high, e.g.
#' `num_runs = 100`.
#'
#' The larger `perc_test` is, the more the test set turns into a matrix
#' completion problem, rather than the desired matrix decomposition problem. To
#' better resemble the actual problem PCP will be faced with come inference
#' time, `perc_test` should therefore be kept relatively low.
#'
#' Choosing a reasonable value for `num_runs` is dependent on the need to keep
#' `perc_test` relatively low. Ideally, a large enough `num_runs` is used so
#' that many (if not all) of the entries in `D` are likely to eventually be
#' tested. Note that since test set entries are chosen randomly for all runs `1`
#' through `num_runs`, in the pathologically worst case scenario, the same exact
#' test set could be drawn each time. In the best case scenario, a different
#' test set is obtained each run, providing balanced coverage of `D`. Viewed
#' another way, the smaller `num_runs` is, the more the results are susceptible
#' to overfitting to the relatively few selected test sets.
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
#'   `max_iter` argument to reach convergence. [rrmc()] does not report
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
#' @param pcp_fn The PCP function to use when grid searching. Must be either
#'   `rrmc` or `root_pcp` (passed without the soft brackets).
#' @param grid A `data.frame` of dimension `j` by `k` containing the `j`-many
#'   unique settings of `k`-many parameters to try.
#'   **NOTE: The columns of `grid` should be
#'   named after the required parameters in the function header of `pcp_fn`.**
#'   For example, if `pcp_fn = root_pcp` and you want to search through `lambda`
#'   and `mu`, then `names(grid)` must be set to `c("lambda", "mu")`. If instead
#'   you want to keep e.g. `lambda` fixed and search through only `mu`, you can
#'   either have a `grid` with only one column, `mu`, and pass `lambda` as a
#'   constant via `...`, or you can have `names(grid)` set to
#'   `c("lambda", "mu")` where `lambda` is constant. The same logic applies for
#'   `pcp_fn = rrmc` and `eta` and `r`.
#' @param ... Any parameters required by `pcp_fn` that should be kept constant
#'   throughout the grid search, or those parameters that cannot be stored in
#'   `grid` (e.g. the `LOD` parameter). A parameter should not be passed with
#'   `...` if it is already a column in `grid`, as that behavior is ambiguous.
#' @param parallel_strategy (Optional) The parallelization strategy used when
#'   conducting the grid search (to be passed on to the [future::plan()]
#'   function). Must be one of: `"sequential"`, `"multisession"`, `"multicore"`
#'   or `"cluster"`. By default, `parallel_strategy = "sequential"`, which
#'   runs the grid search in serial and the `num_workers` argument is ignored.
#'   The option `parallel_strategy = "multisession"` parallelizes the search
#'   via sockets in separate R _sessions_. The option
#'   `parallel_strategy = "multicore"` is not supported on Windows
#'   machines, nor in .Rmd files (must be run in a .R script) but parallelizes
#'   the search much faster than `"multisession"` since it runs separate
#'   _forked_ R processes. The option `parallel_strategy = "cluster"`
#'   parallelizes using separate R sessions running typically on one or more
#'   machines. Support for other parallel strategies will be added in a future
#'   release of `pcpr`. **It is recommended to use
#'   `parallel_strategy = "multicore"` or `"multisession"` when possible.**
#' @param num_workers (Optional) An integer specifying the number of workers to
#'   use when parallelizing the grid search, to be passed on to
#'   [future::plan()]. By default, `num_workers = 1`. When possible, it is
#'   recommended to use `num_workers = parallel::detectCores(logical = F)`,
#'   which computes the number of physical CPUs available on the machine
#'   (see [parallel::detectCores()]). `num_workers` is ignored
#'   when `parallel_strategy = "sequential"`, and must be `> 1` otherwise.
#' @param perc_test (Optional) The fraction of entries of `D` that will be
#'   randomly corrupted as `NA` missing values (the test set). Can be anthing in
#'   the range `[0, 1)`. By default, `perc_test = 0.05`. See **Best practices**
#'   section for more details.
#' @param num_runs (Optional) The number of times to test a given parameter
#'   setting. By default, `num_runs = 100`. See **Best practices** section for
#'   more details.
#' @param return_all_tests (Optional) A logical indicating if you would like the
#'   output from all the calls made to `pcp_fn` over the course of the grid
#'   search to be returned to you in list format. If set to `FALSE`, then only
#'   statistics on the parameters tested will be returned. If set to `TRUE` then
#'   every `L`, and `S` matrix recovered during the grid search will be returned
#'   in the lists `L_mats` and `S_mats`, every test set matrix will be returned
#'   in the list `test_mats`, the original input matrix will be returned as
#'   `original_mat`, and the parameters passed in to `...` will be returned in
#'   the `constant_params` list. **By default, `return_all_tests = FALSE`,
#'   which is highly recommended. Setting `return_all_tests = TRUE` can consume
#'   a massive amount of memory depending on the size of `grid`, the input
#'   matrix `D`, and the value for `num_runs`.**
#' @param verbose (Optional) A logical indicating if you would like verbose
#'   output displayed or not. By default, `verbose = TRUE`. To obtain
#'   progress bar updates, the user must wrap the `grid_search_cv()` call
#'   with a call to [progressr::with_progress()]. The progress bar does _not_
#'   depend on the value passed for `verbose`.
#'
#' @returns A list containing:
#' * `all_stats`: A `data.frame` containing the statistics of every run
#'   comprising the grid search. These statistics include the parameter
#'   settings for the run, along with the `run_num` (used as the seed
#'   for the corruption step, step 1 in the grid search procedure),
#'   the relative error for the run `rel_err`, the rank of the recovered L
#'   matrix `L_rank`, the sparsity of the recovered S matrix `S_sparsity`,
#'   the number of `iterations` PCP took to reach convergence (for [root_pcp()]
#'   only), and the error status `run_error` of the PCP run (`NA` if no error,
#'   otherwise a character string).
#' * `summary_stats`: A `data.frame` containing a summary of the information in
#'   `all_stats`. Summary made by column-wise averaging the results in
#'   `all_stats`.
#' * `metadata`: A character string containing the metadata associated with the
#'   grid search instance.
#'
#' If `return_all_tests = TRUE` then the following are also returned as part
#' of the list:
#' * `L_mats`: A list containing all the `L` matrices returned from PCP
#'   throughout the grid search. Therefore, `length(L_mats) == nrow(all_stats)`.
#'   Row `i` in `all_stats` corresponds to `L_mats[[i]]`.
#' * `S_mats`: A list containing all the S matrices returned from PCP throughout
#'   the grid search. Therefore, `length(S_mats) == nrow(all_stats)`. Row `i` in
#'   `all_stats` corresponds to `S_mats[[i]]`.
#' * `test_mats`: A list of `length(num_runs)` containing all the corrupted test
#'   mats (and their masks) used throughout the grid search. Note:
#'   `all_stats$run[i]` corresponds to `test_mats[[i]]`.
#' * `original_mat`: The original data matrix `D`.
#' * `constant_params`: A copy of the constant parameters that were originally
#'   passed to the grid search (for record keeping).
#' @seealso [sim_na()], [sparsity()], [matrix_rank()], [get_pcp_defaults()]
#' @examples
#' #### -------Simple simulated PCP problem-------####
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
#' gs$summary_stats
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
#' @importFrom rlang .data
grid_search_cv <- function(
    D,
    pcp_fn,
    grid,
    ...,
    parallel_strategy = "sequential",
    num_workers = 1,
    perc_test = 0.05,
    num_runs = 100,
    return_all_tests = FALSE,
    verbose = TRUE) {
  # 0. Error handling:
  constant_params <- list(...)
  constant_params_to_report <- list(...)
  repeated_vars <- intersect(names(constant_params), colnames(grid))
  if (length(repeated_vars) > 0) stop(paste0('Arguments passed to "..." and "grid" are in conflict with one another. The following variables appear in both arguments and are therefore ambiguous: ', paste(repeated_vars, collapse = ", ")))
  if ("lambda" %in% names(constant_params)) {
    grid$lambda <- constant_params$lambda
    constant_params[["lambda"]] <- NULL
  }
  if ("mu" %in% names(constant_params)) {
    grid$mu <- constant_params$mu
    constant_params[["mu"]] <- NULL
  }
  if ("eta" %in% names(constant_params)) {
    grid$eta <- constant_params$eta
    constant_params[["eta"]] <- NULL
  }
  if ("r" %in% names(constant_params)) {
    grid$r <- constant_params$r
    constant_params[["r"]] <- NULL
  }
  checkmate::assert_matrix(D)
  checkmate::assert_choice(parallel_strategy, choices = c("sequential", "multisession", "multicore", "cluster"))
  if (parallel_strategy != "sequential" && num_workers == 1) stop('Arguments "parallel_strategy" and "num_workers" are in conflict with one another. If you want to use 1 core then "parallel_strategy" should be set to "sequential". If you want to run the grid search in parallel then num_workers should be > 1.')
  if (parallel_strategy == "multicore" && Sys.info()["sysname"] == "Windows") stop('Argument "parallel_strategy" cannot be set to "multicore" on a Windows machine. Please use "multisession", "cluster", or "sequential" instead.')
  checkmate::qassert(num_workers, rules = "X1[1,)")
  checkmate::qassert(perc_test, rules = "N1[0,1)")
  checkmate::qassert(num_runs, rules = "X1[1,)")
  if (perc_test == 0 && num_runs > 1) stop('Arguments "perc_test" and "num_runs" are in conflict with one another. Argument "num_runs" cannot be greater than 1 if "perc_test" is 0.')
  checkmate::qassert(return_all_tests, rules = "B1")
  checkmate::qassert(verbose, rules = "B1")
  # 1. Setting up the search:
  if (verbose) cat("\nInitializing grid search...")
  # 1a. Parsing the grid that was passed:
  metrics <- c("rel_err", "L_rank", "S_sparsity", "iterations", "run_error", "run_error_perc")
  for (metric in metrics) {
    if (!metric %in% names(grid)) {
      grid[, metric] <- NA
    }
  }
  # 1b. Setting up the params to search through:
  param_names <- grid %>%
    dplyr::select(!tidyselect::all_of(metrics)) %>%
    colnames()
  points_to_eval <- which(is.na(grid$rel_err))
  params <- data.frame(grid[points_to_eval, param_names], run_num = rep(1:num_runs, each = length(points_to_eval)), row.names = NULL)
  colnames(params) <- c(param_names, "run_num")
  # 1c. Creating the test sets:
  test_mats <- purrr::map(1:num_runs, ~ sim_na(D, perc = perc_test, seed = .x))
  # 1d. Parallel programming setup:
  start_time <- Sys.time()
  if (num_workers == 1) {
    future::plan(parallel_strategy)
    if (verbose) cat(paste0("\nBeginning sequential grid search...\nStart time: ", start_time, "\n"))
  } else {
    future::plan(parallel_strategy, workers = num_workers)
    if (verbose) cat(paste0("\nBeginning parallel grid search using ", num_workers, " cores and a ", parallel_strategy, " strategy...\nStart time: ", start_time, "\n"))
  }
  # 1e. Progress bar setup:
  p <- progressr::progressor(steps = nrow(params))
  # 2. Conducting the search:
  if (!return_all_tests) {
    pcp_evals <- furrr::future_map_dfr(1:nrow(params), .f = function(i) {
      p()
      tryCatch(expr = {
        param_setting <- as.data.frame(params[i, param_names])
        colnames(param_setting) <- param_names
        pcp_mat <- do.call(pcp_fn, c(constant_params, as.list(param_setting), list(D = test_mats[[params$run_num[i]]]$D_tilde)))
        eval_params(settings = params[i, ], test_mat = D, pcp_model = pcp_mat, test_mask = test_mats[[params$run_num[i]]]$tilde_mask)
      }, error = function(e) {
        cbind(params[i, ], data.frame(rel_err = NA, L_rank = NA, S_sparsity = NA, iterations = NA, run_error = paste0(e)))
      })
    })
    future::plan("sequential")
  } else {
    pcp_mats <- furrr::future_map(1:nrow(params), .f = function(i) {
      p()
      param_setting <- as.data.frame(params[i, param_names])
      colnames(param_setting) <- param_names
      tryCatch(expr = {
        do.call(pcp_fn, c(constant_params, as.list(param_setting), list(D = test_mats[[params$run_num[i]]]$D_tilde)))
      }, error = function(e) {
        list(message = paste0(e), settings = as.list(param_setting), test_mat = test_mats[[params$run_num[i]]]$D_tilde)
      })
    })
    future::plan("sequential")
    pcp_evals <- purrr::map_dfr(1:nrow(params), .f = function(i) {
      eval_params(settings = params[i, ], test_mat = D, pcp_model = pcp_mats[[i]], test_mask = test_mats[[params$run_num[i]]]$tilde_mask)
    })
  }
  end_time <- Sys.time()
  if (verbose) cat(paste0("\nGrid search completed at time: ", end_time))
  if (perc_test == 0) pcp_evals$rel_err <- NA
  # 3. Summarizing the results from the search:
  if (num_runs > 1) {
    evals_summary <- pcp_evals %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(param_names))) %>%
      dplyr::summarise(rel_err = mean(.data[["rel_err"]], na.rm = T), L_rank = mean(.data[["L_rank"]], na.rm = T), S_sparsity = mean(.data[["S_sparsity"]], na.rm = T), iterations = mean(.data[["iterations"]], na.rm = T), run_error_perc = paste0(round(100 * sum(!is.na(.data[["run_error"]])) / dplyr::n(), 2), "%"), .groups = "drop") %>%
      dplyr::arrange(.data[["rel_err"]])
  } else {
    evals_summary <- pcp_evals %>%
      dplyr::select(!.data[["run_num"]]) %>%
      dplyr::arrange(.data[["rel_err"]])
  }
  if (verbose) cat("\nMetrics calculations complete.")
  # 4. package results:
  if (!return_all_tests) {
    results <- list(all_stats = pcp_evals, summary_stats = evals_summary)
  } else {
    if ("L_list" %in% names(pcp_mats[[1]])) {
      Ls <- foreach::`%do%`(
        foreach::foreach(k = 1:length(pcp_mats), .combine = c),
        pcp_mats$k$L_list
      )
      Ss <- foreach::`%do%`(
        foreach::foreach(k = 1:length(pcp_mats), .combine = c),
        pcp_mats$k$S_list
      )
    } else {
      Ls <- foreach::`%do%`(
        foreach::foreach(k = 1:length(pcp_mats), .combine = c),
        pcp_mats$k$L
      )
      Ss <- foreach::`%do%`(
        foreach::foreach(k = 1:length(pcp_mats), .combine = c),
        pcp_mats$k$S
      )
    }
    results <- list(all_stats = pcp_evals, summary_stats = evals_summary, L_mats = Ls, S_mats = Ss, test_mats = test_mats, original_mat = D, constant_params = constant_params_to_report)
  }
  # 5. metadata to return to user:
  params_summary <- c()
  for (p in param_names) params_summary <- c(params_summary, paste0("\t\t", p, ": {", paste(sort(unique(params[[p]])), collapse = ", "), "}"))
  metadata <- c(
    paste(rep("-", 80), collapse = ""),
    paste0("Date & Time: ", Sys.time()),
    paste0("Description: Metadata string for pcpr's grid_search_cv()\n"),
    paste(rep("-", 80), collapse = ""),
    paste0("\tGrid search start time: ", start_time),
    paste0("\tGrid search end time: ", end_time),
    paste0("\tPCP function used: ", deparse(pcp_fn)[1]),
    paste0("\tTotal number of settings searched through: ", nrow(params)),
    paste0("\tParameters searched through: "),
    params_summary,
    paste0("\nGrid search settings:"),
    paste0("\tParallelization strategy used: ", parallel_strategy),
    paste0("\tCores used: ", num_workers),
    paste0("\tPercent of matrix corrupted for test set: ", round(perc_test * 100, 3), "%"),
    paste0("\tRuns per parameter setting: ", num_runs),
    paste0("\tReturn all matrices from all tests?: ", return_all_tests),
    paste0("\tConstant parameters passed to PCP: ", paste(names(constant_params_to_report), collapse = ", "), "\n"),
    paste(rep("-", 80), collapse = "")
  )
  if (verbose) cat("\nGrid search completed!")
  results$metadata <- metadata
  # 6. Return results:
  results
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
      rel_err = norm((test_mat - pcp_model$L_list[[.x]]) * test_mask, type = "F") / norm(test_mat * test_mask, type = "F"),
      L_rank = matrix_rank(pcp_model$L_list[[.x]], 1e-04),
      S_sparsity = sparsity(pcp_model$S_list[[.x]]),
      iterations = NA, run_error = NA
    ))
  } else {
    test_mat[is.na(test_mat)] <- 0
    stats <- settings
    stats$rel_err <- norm((test_mat - pcp_model$L) * test_mask, type = "F") / norm(test_mat * test_mask, type = "F")
    stats$L_rank <- matrix_rank(pcp_model$L, 1e-04)
    stats$S_sparsity <- sparsity(pcp_model$S)
    stats$iterations <- pcp_model$num_iter
    stats$run_error <- NA
  }
  stats
}
