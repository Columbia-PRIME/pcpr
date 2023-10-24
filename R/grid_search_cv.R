#' Cross-validated grid search of the parameters for Principle Component Pursuit (PCP).
#'
#' \code{grid_search_cv} conducts a cross-validated grid search of the parameters for a given data matrix \code{mat}, PCP function \code{pcp_func}, and settings of parameters to search through \code{grid}. See the \strong{Methods} section below for more details.
#'
#' @param mat The data matrix to conduct the grid search on.
#' @param pcp_func The PCP function to use when grid searching. \emph{Note: the PCP function passed must be able to handle missing \code{NA} values.} For example: \code{root_pcp_na}.
#' @param grid A dataframe with dimension N x P containing the N-many settings of P-many parameters to try. \emph{The columns of \code{grid} should be named exactly as they are in the function
#' header of} \code{pcp_func}. For example, if \code{pcp_func = root_pcp_noncvx_na}, then the columns of \code{grid} should be named "lambda", "mu", and "r"
#' An optional additional column named "rel_err" can be included that contains the mean relative error recorded by that row's parameter setting, with those rows (settings) that have not been tried left as \code{NA}.
#' In this way, you can perform a grid search in which you already know the relative errors of some parameter settings, but would like to expand your knowledge of the unexplored parts of the grid further.
#' @param scale_func (Optional) The function used to scale the input \code{mat} by column. By default, \code{scale_func = NULL}, and no scaling is to be done at all.
#' @param parallel_approach (Optional) The computational approach used when conducting the gridsearch (to be passed on to the \code{future} package's \code{\link[future]{plan}} function). Must be one of: \code{"sequential", "multisession", "multicore"}. By default, \code{parallel_approach = "multisession"}, which does parallelization via sockets (in separate R sessions) and works on any operating system. If \code{parallel_approach = "sequential"} then the search will be conducted in serial. The option \code{parallel_approach = "multicore"} is not supported on Windows machines, nor in RStudio (must be run from the command line) but is faster than the \code{"multisession"} approach since it runs separate forked R processes.
#' @param cores (Optional) The number of cores to use when parallelizing the grid search. By default, \code{cores = parallel::detectCores(logical = F)}, which is the number of physical CPUs available on the machine.
#' @param perc_test (Optional) The fraction of entries of \code{mat} that will be randomly imputed as \code{NA} missing values (the test set). Can be anthing in the range \code{[0, 1)}. By default, \code{perc_test = 0.15}.
#' @param runs (Optional) The number of times to test a given parameter setting. By default, \code{runs = 1}.
#' @param conserve_memory (Optional) A logical indicating if you only care about the actual statistics of the gridsearch and would therefore like to conserve memory when running the gridsearch. If set to \code{TRUE}, then only statistics on the parameters tested will be returned. By default, \code{conserve_memory = FALSE}, in which case additional objects saving the outputs of all runs of \code{pcp_func} will also be returned.
#' @param verbose (Optional) A logical indicating if you would like verbose output displayed or not. By default, \code{verbose = TRUE}.
#' @param save_as (Optional) A character containing the root of the file path used to save the output to. Importantly, this should not end in any file extension, since this character will be used to save both the resulting \code{[save_as].rds} and {[save_as]_README.txt} files. By default, \code{save_as = NULL}, in which case the gridsearch is not saved to any file.
#' @param ... \emph{Any parameters required by \code{pcp_func} that could not be specified in \code{grid}. Importantly, these parameters are therefore kept constant (not involved in the grid search). The best example is the \code{LOD} parameter for those PCP functions that require the \code{LOD} argument.}
#'
#' @section Methods:
#' Each hyperparameter setting is cross-validated by:
#' \enumerate{
#'   \item Randomly corrupting \code{perc_test} percent of the entries in \code{mat} as missing (i.e. \code{NA} values), yielding \code{corrupted_mat}.
#'   Done via the \code{\link{corrupt_mat_randomly}} function.
#'   \item Running the PCP function (\code{pcp_func}) on \code{corrupted_mat}, giving \code{L_hat} and \code{S_hat}.
#'   \item Recording the relative recovery errors of \code{L_hat} compared with the input data matrix \code{mat} for only those values that were imputed as missing during the corruption step. Ie. \code{||P_OmegaCompliment(mat - L_hat)||_F / ||P_OmegaCompliment(mat)||_F}.
#'   \item Repeating steps 1-3 for a total of \code{runs} many times.
#'   \item Reporting the mean of the \code{runs}-many runs for each parameter setting.
#' }
#'
#' @return A list containing the following:
#' \describe{
#'    \item{\code{all_stats}}{A \code{data.frame} containing the statistics of every run comprising the grid search. These statistics include the parameter settings for the run, along with the \code{run} number (used as the seed in the corruption step outlined in step 1 of the \strong{Methods} section), the relative error for the run \code{rel_err}, the rank of the recovered L matrix \code{L_rank}, the sparsity of the recovered S matrix \code{S_sparsity}, the number of \code{iterations} PCP took to reach convergence, and the error status \code{run_error} of the PCP run (\code{NA} if no error, otherwise a character).}
#'    \item{\code{summary_stats}}{A \code{data.frame} containing a summary of the information in \code{all_stats}. Made to easily pass on to \code{\link{print_gs}}.}
#'    \item{\code{L_mats}}{A list containing all the L matrices returned from PCP throughout the gridsearch. Therefore, \code{length(L_mats) == nrow(all_stats)}. Row i in \code{all_stats} corresponds to \code{L_mats[[i]]}. Only returned when \code{conserve_memory = FALSE}.}
#'    \item{\code{S_mats}}{A list containing all the S matrices returned from PCP throughout the gridsearch. Therefore, \code{length(S_mats) == nrow(all_stats)}. Row i in \code{all_stats} corresponds to \code{S_mats[[i]]}. Only returned when \code{conserve_memory = FALSE}.}
#'    \item{\code{test_mats}}{A list of \code{length(runs)} containing all the corrupted test mats (and their masks) used throughout the gridsearch. Note: \code{all_stats$run[i]} corresponds to \code{test_mats[[i]]}. Only returned when \code{conserve_memory = FALSE}.}
#'    \item{\code{original_mat}}{The original data matrix \code{mat} after it was column scaled by \code{scale_func}. Only returned when \code{conserve_memory = FALSE}.}
#'    \item{\code{constant_params}}{A copy of the constant parameters that were originally passed to the gridsearch (for record keeping).}
#' }
#' @seealso Older versions of PCP's gridsearch (not recommended): \code{\link{grid_search_cv}}, \code{\link{random_search_cv}}, \code{\link{bayes_search_cv}}, and \code{\link{print_gs}}
#' @examples
#'
#' library(pcpr) # since we will be passing \code{grid_search_cv} a PCP function
#'
#' # simulate a data matrix:
#'
#' n <- 50
#' p <- 10
#' data <- sim_data(sim_seed = 1, nrow = n, ncol = p, rank = 3, sigma=0, add_sparse = FALSE)
#' mat <- data$M
#'
#' # pick parameter settings of lambda and mu to try:
#'
#' lambdas <- c(1/sqrt(n), 1.25/sqrt(n), 1.5/sqrt(n))
#' mus <- c(sqrt(p/2), sqrt(p/1.5), sqrt(p/1.25))
#' param_grid <- expand.grid(lambda = lambdas, mu = mus)
#'
#' # run the grid search:
#'
#' search_results <- grid_search_cv(mat, pcp_func = root_pcp_na, grid_df = param_grid, cores = 4, perc_b = 0.2, runs = 20, verbose = TRUE, save_as = NULL)
#'
#' # visualize the output:
#'
#' print_gs2(search_results$summary_stats)
#' @export
#' @importFrom magrittr %>%
grid_search_cv <- function(
  mat,
  pcp_func,
  grid,
  scale_func = NULL,
  parallel_approach = "multisession",
  cores = parallel::detectCores(logical = F),
  perc_test = 0.05,
  runs = 100,
  conserve_memory = FALSE,
  verbose = TRUE,
  save_as = NULL,
  ...
)
{

  # 0. Error handling:
  constant_params <- list(...)
  repeated_vars <- intersect(names(constant_params), colnames(grid))
  if (length(repeated_vars) > 0) stop(paste0('Arguments passed to "..." and "grid" are in conflict with one another. The following variables appear in both arguments and are therefore ambiguous: ', paste(repeated_vars, collapse = ", ")))
  if (!is.matrix(mat)) stop('Invalid value passed for argument "mat". Argument "mat" must be a matrix. See documentation with "?grid_search_cv".')
  if (!(parallel_approach %in% c("sequential", "multisession", "multicore"))) stop('Invalid value passed for argument "parallel_approach". Must be one of: {"sequential", "multisession", "multicore"}. See documentation with "?grid_search_cv".')
  if (parallel_approach == "sequential" && cores != 1) stop('Arguments "parallel_approach" and "cores" are in conflict with one another. If you want the "parallel_approach" to be "sequential" then cores should be set to 1. If you want to use more than 1 core then "parallel_approach" should be set as either "multisession" or "multicore".')
  if (parallel_approach != "sequential" && cores == 1) stop('Arguments "parallel_approach" and "cores" are in conflict with one another. If you want to use 1 core then "parallel_approach" should be set to "sequential". If you want to run the gridsearch in parallel then cores should be > 1.')
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
  if (!is.null(scale_func)) mat <- apply(mat, 2, scale_func)

  # 1e. Creating the test sets:
  test_mats <- purrr::map(1:runs, ~ corrupt_mat_randomly(mat, seed = .x, perc = perc_test))

  # 1f. Parallel programming setup:
  start_time <- Sys.time()
  if (cores == 1) {
    future::plan(parallel_approach)
    if (verbose) cat(paste0("\nBeginning sequential gridsearch...\nStart time: ", start_time, "\n"))
  } else {
    future::plan(parallel_approach, workers = cores)
    if (verbose) cat(paste0("\nBeginning parallel gridsearch using ", cores, " cores and a ", parallel_approach, " approach...\nStart time: ", start_time, "\n"))
  }

  # 1g. Progress bar setup:
  p <- progressr::progressor(steps = nrow(params))
  old_progress_bar <- verbose && parallel_approach == "multicore"

  # 2. Conducting the search:
  if (conserve_memory) {
    pcp_evals <- furrr::future_map_dfr(1:nrow(params), .f = function(i) {
      p()
      tryCatch(expr = {
        param_setting <- as.data.frame(params[i, param_names])
        colnames(param_setting) <- param_names
        pcp_mat <- do.call(pcp_func, c(constant_params, as.list(param_setting), list(D = test_mats[[params$run_num[i]]]$mat_tilde)))
        eval_params(settings = params[i,], test_mat = mat, pcp_mat = pcp_mat, test_mask = test_mats[[params$run_num[i]]]$tilde_mask)
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
        do.call(pcp_func, c(constant_params, as.list(param_setting), list(D = test_mats[[params$run_num[i]]]$mat_tilde)))
      }, error = function(e) {
        list(message = paste0(e), settings = as.list(param_setting), test_mat = test_mats[[params$run_num[i]]]$mat_tilde)
      })
    }, .progress = old_progress_bar)

    future::plan("sequential")
    pcp_evals <- purrr::map_dfr(1:nrow(params), .f = function(i) {
      eval_params(settings = params[i,], test_mat = mat, pcp_mat = pcp_mats[[i]], test_mask = test_mats[[params$run_num[i]]]$tilde_mask)
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

    results <- list(all_stats = pcp_evals, summary_stats = evals_summary, L_mats = Ls, S_mats = Ss, test_mats = test_mats, original_mat = mat, constant_params = constant_params)

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
      paste0("\nGridsearch type: vanilla"),
      paste0("\tGridsearch start time: ", start_time),
      paste0("\tGridsearch end time: ", end_time),
      paste0("\tPCP function used: ", deparse(pcp_func)[1]),
      paste0("\tTotal number of settings searched through: ", nrow(params)),
      paste0("\tParameters searched through: "),
      params_summary,
      paste0("\nGridsearch settings:"),
      paste0("\tScale function used: ", deparse(scale_func)[1]),
      paste0("\tParallelization approach used: ", parallel_approach),
      paste0("\tCores used: ", cores),
      paste0("\tPercent of matrix corrupted for test set: ", round(perc_test*100, 3), "%"),
      paste0("\tRuns per parameter setting: ", runs),
      paste0("\tMemory conserved?: ", conserve_memory),
      paste0("\tConstant parameters passed to PCP: ", paste(names(constant_params), collapse = ", "), "\n"),
      paste(rep('-', 80), collapse ="")
    ), README)
    close(README)
    if (verbose) cat("\nSave completed. Exiting gridsearch.")
  }

  # 6. Return results:
  results

}

#' @keywords internal
corrupt_mat_randomly <- function(mat, seed, perc) {
  nvals_to_corrupt <- floor(perc*prod(dim(mat)))
  mat_vec <- as.vector(mat)
  mask <- rep(0, length(mat_vec))

  pool <- which(mat_vec >= 0)
  if (length(pool) == 0) stop('There is nothing in the input matrix "mat" that can be corrupted as missing.')
  if (nvals_to_corrupt > length(pool)) {
    corrupted <- pool
  } else {
    set.seed(seed)
    corrupted <- sample(pool, nvals_to_corrupt, replace = FALSE)
  }

  mask[corrupted] <- 1
  mat_vec[corrupted] <- NA

  rows <- nrow(mat)
  cols <- ncol(mat)

  ret_mat <- matrix(mat_vec, nrow = rows, ncol = cols)
  ret_mask <- matrix(mask, nrow = rows, ncol = cols)

  list(mat_tilde = ret_mat, tilde_mask = ret_mask)
}

#' @keywords internal
eval_params <- function(settings, test_mat, pcp_mat, test_mask) {

  if ("message" == names(pcp_mat)[1]) {

    stats <- cbind(settings, data.frame(rel_err = NA, L_rank = NA, S_sparsity = NA, iterations = NA, run_error = pcp_mat$message))

  } else if ("L_list" == names(pcp_mat)[3]) {

    if ("r" %in% names(settings)) settings <- settings[setdiff(names(settings), "r")]
    test_mat[is.na(test_mat)] <- 0
    stats <- purrr::imap_dfr(.x = 1:length(pcp_mat$L_list), ~ data.frame(
      r = .y, settings,
      rel_err = norm((test_mat - pcp_mat$L_list[[.x]])*test_mask, "F") / norm(test_mat*test_mask, "F"),
      L_rank = Matrix::rankMatrix(pcp_mat$L_list[[.x]], tol = 1e-04),
      S_sparsity = sparsity(pcp_mat$S_list[[.x]], c = 1e-04),
      iterations = NA, run_error = NA
    ))

  } else {

    test_mat[is.na(test_mat)] <- 0
    stats <- settings
    stats$rel_err <- norm((test_mat - pcp_mat$L)*test_mask, "F") / norm(test_mat*test_mask, "F")
    stats$L_rank <- Matrix::rankMatrix(pcp_mat$L, tol = 1e-04)
    stats$S_sparsity <- sparsity(pcp_mat$S, c = 1e-04)
    stats$iterations <- pcp_mat$final_iter
    stats$run_error <- NA

  }

  stats

}

#' Calculate the sparsity of a given matrix.
#'
#' \code{sparsity} computes the sparsity of a given matrix, \code{S}, using a given tolerance, \code{c}.
#'
#' @param S The input matrix whose sparsity is in question.
#' @param c (Optional) a numeric threshold / tolerance used to determine if an entry in \code{S}
#' is effectively 0. If the absolute value of an entry is below \code{c} then it is effectively 0.
#' By default, \code{c = 1e-05}.
#'
#' @return The sparsity of \code{S}.
#' @examples
#' sparsity(matrix(0:8, 3, 3))
#' @export
sparsity <- function(S, c = 1e-05) {
  S[abs(S) < c] <- 0
  100 - (100 * sum((S != 0) * 1) / prod(dim(S)))
}
