#!/usr/bin/env Rscript

# Compare the current local hOUwie implementation with phytools::fitmultiOU.
#
# Usage:
#   Rscript tools/benchmarks/houwie-fitmultiou/run-comparison.R smoke
#   Rscript tools/benchmarks/houwie-fitmultiou/run-comparison.R pilot
#   Rscript tools/benchmarks/houwie-fitmultiou/run-comparison.R scaling
#   Rscript tools/benchmarks/houwie-fitmultiou/run-comparison.R article
#   Rscript tools/benchmarks/houwie-fitmultiou/run-comparison.R study
#
# fitmultiOU is currently available in the development version of phytools. If
# that version is installed in a separate library, point this script to it with:
#   PHYTOOLS_LIB=/path/to/library Rscript ... pilot
#
# Results are written after every task, so an interrupted run can be resumed.
# Set HOUWIE_BENCH_OVERWRITE=1 to replace results for the selected profile.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
benchmark_dir <- dirname(script_file)
repo_root <- normalizePath(file.path(benchmark_dir, "..", "..", ".."),
                           mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
profile_name <- if(length(args)) args[1] else "pilot"

profiles <- list(
  smoke = list(
    tips = 30L,
    data_seeds = 8L,
    fit_repetitions = 1L,
    evaluation_repetitions = 1L,
    houwie_nsim = 8L,
    fitmultiou_levels = 15L,
    houwie_eval_nsim = c(4L, 8L),
    fitmultiou_eval_levels = c(10L, 15L),
    houwie_maxeval = 25L,
    fitmultiou_maxit = 25L,
    houwie_starts = 1L,
    houwie_use_matched_init = TRUE,
    include_parallel_houwie = FALSE,
    default_cores = 1L
  ),
  pilot = list(
    tips = 75L,
    data_seeds = 8L,
    fit_repetitions = 3L,
    evaluation_repetitions = 5L,
    houwie_nsim = 25L,
    fitmultiou_levels = 40L,
    houwie_eval_nsim = c(10L, 25L, 50L),
    fitmultiou_eval_levels = c(20L, 40L, 80L),
    houwie_maxeval = 200L,
    fitmultiou_maxit = 200L,
    houwie_starts = 1L,
    houwie_use_matched_init = TRUE,
    include_parallel_houwie = FALSE,
    default_cores = 1L
  ),
  scaling = list(
    tips = c(50L, 100L, 300L),
    data_seeds = 8L,
    fit_repetitions = 0L,
    evaluation_repetitions = 3L,
    houwie_nsim = 100L,
    fitmultiou_levels = 100L,
    houwie_eval_nsim = c(25L, 100L),
    fitmultiou_eval_levels = c(50L, 100L),
    houwie_maxeval = 1000L,
    fitmultiou_maxit = 1000L,
    houwie_starts = 1L,
    houwie_use_matched_init = TRUE,
    include_parallel_houwie = FALSE,
    default_cores = 1L
  ),
  article = list(
    tips = 300L,
    data_seeds = 8L,
    fit_repetitions = 1L,
    evaluation_repetitions = 0L,
    houwie_nsim = 100L,
    fitmultiou_levels = 100L,
    houwie_eval_nsim = integer(0),
    fitmultiou_eval_levels = integer(0),
    houwie_maxeval = 1000L,
    fitmultiou_maxit = 2000L,
    houwie_starts = 10L,
    houwie_use_matched_init = FALSE,
    include_parallel_houwie = TRUE,
    default_cores = 10L
  ),
  study = list(
    tips = c(50L, 100L, 300L),
    data_seeds = c(8L, 18L, 28L),
    fit_repetitions = 3L,
    evaluation_repetitions = 3L,
    houwie_nsim = 100L,
    fitmultiou_levels = 100L,
    houwie_eval_nsim = c(25L, 50L, 100L, 200L),
    fitmultiou_eval_levels = c(25L, 50L, 100L, 200L),
    houwie_maxeval = 1000L,
    fitmultiou_maxit = 1000L,
    houwie_starts = 1L,
    houwie_use_matched_init = TRUE,
    include_parallel_houwie = FALSE,
    default_cores = 1L
  )
)

if(!profile_name %in% names(profiles)){
  stop(
    "Unknown profile: ", profile_name,
    ". Choose smoke, pilot, scaling, article, or study.",
    call. = FALSE
  )
}
profile <- profiles[[profile_name]]

requested_cores <- suppressWarnings(as.integer(
  Sys.getenv("HOUWIE_BENCH_CORES", profile$default_cores)
))
if(is.na(requested_cores) || requested_cores < 1L){
  stop("HOUWIE_BENCH_CORES must be a positive integer.", call. = FALSE)
}
physical_cores <- parallel::detectCores(logical = FALSE)
available_cores <- parallel::detectCores(logical = TRUE)
if(is.na(available_cores)) available_cores <- 1L
if(is.na(physical_cores)) physical_cores <- available_cores
cores <- min(requested_cores, available_cores)

phytools_lib <- Sys.getenv("PHYTOOLS_LIB")
if(nzchar(phytools_lib)){
  if(!dir.exists(phytools_lib)){
    stop("PHYTOOLS_LIB does not exist: ", phytools_lib, call. = FALSE)
  }
  .libPaths(c(normalizePath(phytools_lib), .libPaths()))
}

if(!requireNamespace("phytools", quietly = TRUE)){
  stop("phytools is required.", call. = FALSE)
}
if(!"fitmultiOU" %in% getNamespaceExports("phytools")){
  stop(
    "The installed phytools version does not export fitmultiOU. ",
    "Install the development version described in Liam Revell's 31 July 2026 post, ",
    "or set PHYTOOLS_LIB to a library containing it.",
    call. = FALSE
  )
}
if(!requireNamespace("pkgload", quietly = TRUE)){
  stop("pkgload is required to benchmark the current local OUwie source.",
       call. = FALSE)
}

pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                  helpers = FALSE)

results_dir <- file.path(benchmark_dir, "results", profile_name)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

overwrite <- identical(Sys.getenv("HOUWIE_BENCH_OVERWRITE"), "1")
fits_path <- file.path(results_dir, "fits.csv")
evaluations_path <- file.path(results_dir, "evaluations.csv")
if(overwrite){
  unlink(c(fits_path, evaluations_path))
}

read_results <- function(path){
  if(file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
}

append_result <- function(row, path){
  old <- read_results(path)
  combined <- if(is.null(old)) row else rbind(old, row)
  write.csv(combined, path, row.names = FALSE)
  invisible(combined)
}

task_finished <- function(task_id, path){
  old <- read_results(path)
  !is.null(old) && task_id %in% old$task_id
}

simulate_fixture <- function(n_tips, seed){
  set.seed(seed)
  phy <- phytools::pbtree(n = n_tips, scale = 10)
  states <- letters[1:2]
  q <- 0.3
  Q <- matrix(q, 2, 2, dimnames = list(states, states))
  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)
  root_state <- sample(states, 1)
  sim_tree <- phytools::sim.history(phy, Q, anc = root_state)
  alpha <- setNames(rep(0.5, 2), states)
  sigma_sq <- setNames(rep(0.2, 2), states)
  theta <- setNames(c(5, 10), states)
  trait <- phytools::multiOU(
    sim_tree,
    alpha,
    sigma_sq,
    theta,
    a0 = theta[root_state]
  )
  discrete <- factor(
    phytools::getStates(sim_tree, "tips"),
    levels = states
  )
  if(any(table(discrete) == 0)){
    stop("Simulation seed ", seed, " did not retain both states at the tips.",
         call. = FALSE)
  }
  initial <- setNames(
    c(
      mean(trait[discrete == states[1]]),
      mean(trait[discrete == states[2]]),
      log(2) / max(phytools::nodeHeights(phy)),
      var(trait) / max(phytools::nodeHeights(phy)),
      phytools::fitMk(phy, discrete, model = "ER")$rates
    ),
    c("theta[a]", "theta[b]", "alpha", "sigsq", "q[1]")
  )
  data <- data.frame(
    species = names(trait),
    state = discrete,
    trait = as.numeric(trait),
    stringsAsFactors = FALSE
  )
  list(
    phy = phy,
    sim_tree = sim_tree,
    trait = trait,
    discrete = discrete,
    data = data,
    initial = initial,
    truth = c(theta_a = 5, theta_b = 10, alpha = 0.5,
              sigma_sq = 0.2, q = 0.3),
    root_state = root_state,
    seed = seed
  )
}

houwie_initial <- function(fixture){
  unname(fixture$initial[c("q[1]", "alpha", "sigsq", "theta[a]", "theta[b]")])
}

extract_houwie <- function(fit){
  state_columns <- paste0("(", seq_len(2), ")")
  c(
    loglik = fit$loglik,
    theta_a = fit$solution.cont["theta", state_columns[1]],
    theta_b = fit$solution.cont["theta", state_columns[2]],
    alpha = fit$solution.cont["alpha", state_columns[1]],
    sigma_sq = fit$solution.cont["sigma2", state_columns[1]],
    q = mean(fit$solution.disc[!is.na(fit$solution.disc)]),
    evaluations = if(!is.null(fit$global_liks_mat))
      sum(fit$global_liks_mat[[1]] != 0) else NA_real_,
    convergence = NA_real_
  )
}

extract_fitmultiou <- function(fit){
  c(
    loglik = as.numeric(fit$logLik),
    theta_a = unname(fit$theta[1]),
    theta_b = unname(fit$theta[2]),
    alpha = unname(fit$alpha),
    sigma_sq = unname(fit$sigsq),
    q = mean(unname(fit$rates)),
    evaluations = unname(fit$opt_results$counts[["function"]]),
    convergence = unname(fit$opt_results$convergence)
  )
}

parameter_errors <- function(estimates, truth){
  theta_rmse <- sqrt(mean(
    (estimates[c("theta_a", "theta_b")] - truth[c("theta_a", "theta_b")])^2
  ))
  c(
    theta_rmse = theta_rmse,
    alpha_log_error = unname(abs(log(estimates["alpha"] / truth["alpha"]))),
    sigma_log_error = unname(abs(log(
      estimates["sigma_sq"] / truth["sigma_sq"]
    ))),
    q_log_error = unname(abs(log(estimates["q"] / truth["q"])))
  )
}

run_fit <- function(method, fixture, run_seed){
  set.seed(run_seed)
  started <- proc.time()[["elapsed"]]
  error_message <- NA_character_
  fit <- tryCatch(
    {
      if(grepl("^hOUwie", method)){
        houwie_cores <- if(method == "hOUwie parallel"){
          min(cores, profile$houwie_starts)
        }else{
          1L
        }
        hOUwie(
          phy = fixture$phy,
          data = fixture$data,
          rate.cat = 1,
          discrete_model = "ER",
          continuous_model = "OUM",
          nSim = profile$houwie_nsim,
          ip = if(profile$houwie_use_matched_init)
            houwie_initial(fixture) else NULL,
          opts = list(
            algorithm = "NLOPT_LN_SBPLX",
            maxeval = profile$houwie_maxeval,
            ftol_rel = .Machine$double.eps^0.25
          ),
          common_random_numbers = TRUE,
          n_starts = profile$houwie_starts,
          ncores = houwie_cores,
          quiet = TRUE,
          sample_tips = FALSE,
          sample_nodes = FALSE,
          adaptive_sampling = FALSE
        )
      }else{
        phytools::fitmultiOU(
          fixture$phy,
          fixture$trait,
          fixture$discrete,
          model = "ER",
          levs = profile$fitmultiou_levels,
          parallel = cores > 1L,
          ncores = cores,
          root = "mle",
          trace = 0,
          maxit = profile$fitmultiou_maxit,
          rand_start = FALSE,
          init = fixture$initial
        )
      }
    },
    error = function(e){
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if(is.null(fit)){
    estimates <- setNames(rep(NA_real_, 8), c(
      "loglik", "theta_a", "theta_b", "alpha", "sigma_sq", "q",
      "evaluations", "convergence"
    ))
    }else if(grepl("^hOUwie", method)){
    estimates <- extract_houwie(fit)
  }else{
    estimates <- extract_fitmultiou(fit)
  }
  errors <- if(all(is.finite(estimates[c(
    "theta_a", "theta_b", "alpha", "sigma_sq", "q"
  )]))) parameter_errors(estimates, fixture$truth) else
    c(theta_rmse = NA_real_, alpha_log_error = NA_real_,
      sigma_log_error = NA_real_, q_log_error = NA_real_)
  c(
    elapsed_seconds = elapsed,
    estimates,
    errors,
    error = error_message
  )
}

run_evaluation <- function(method, resolution, fixture, run_seed){
  set.seed(run_seed)
  started <- proc.time()[["elapsed"]]
  error_message <- NA_character_
  loglik <- tryCatch(
    {
      if(method == "hOUwie"){
        fit <- hOUwie(
          phy = fixture$phy,
          data = fixture$data,
          rate.cat = 1,
          discrete_model = "ER",
          continuous_model = "OUM",
          nSim = resolution,
          p = unname(c(0.3, 0.5, 0.2, 5, 10)),
          common_random_numbers = TRUE,
          quiet = TRUE,
          sample_tips = FALSE,
          sample_nodes = FALSE,
          adaptive_sampling = FALSE
        )
        fit$loglik
      }else{
        initial <- setNames(
          c(5, 10, 0.5, 0.2, 0.3),
          names(fixture$initial)
        )
        fit <- phytools::fitmultiOU(
          fixture$phy,
          fixture$trait,
          fixture$discrete,
          model = "ER",
          levs = resolution,
          parallel = FALSE,
          root = "mle",
          trace = 0,
          maxit = 0,
          rand_start = FALSE,
          init = initial
        )
        as.numeric(fit$logLik)
      }
    },
    error = function(e){
      error_message <<- conditionMessage(e)
      NA_real_
    }
  )
  elapsed <- proc.time()[["elapsed"]] - started
  c(elapsed_seconds = elapsed, loglik = loglik, error = error_message)
}

metadata <- data.frame(
  profile = profile_name,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string,
  ouwie_version = as.character(utils::packageVersion("OUwie")),
  phytools_version = as.character(utils::packageVersion("phytools")),
  phytools_source_sha = Sys.getenv("PHYTOOLS_SOURCE_SHA", NA_character_),
  os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
  physical_cores = physical_cores,
  logical_cores = available_cores,
  benchmark_cores = cores,
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(results_dir, "metadata.csv"), row.names = FALSE)

for(n_tips in profile$tips){
  for(data_seed in profile$data_seeds){
    fixture_path <- file.path(
      results_dir,
      sprintf("fixture-tips%d-seed%d.rds", n_tips, data_seed)
    )
    if(file.exists(fixture_path) && !overwrite){
      fixture <- readRDS(fixture_path)
    }else{
      fixture <- simulate_fixture(n_tips, data_seed)
      saveRDS(fixture, fixture_path)
    }

    for(resolution in profile$houwie_eval_nsim){
      for(evaluation_repetition in seq_len(profile$evaluation_repetitions)){
        task_id <- sprintf("%s-eval-houwie-t%d-d%d-r%d-rep%d", profile_name,
                           n_tips, data_seed, resolution,
                           evaluation_repetition)
        if(!task_finished(task_id, evaluations_path)){
          cat("Running", task_id, "\n")
          run_seed <- data_seed * 10000L + evaluation_repetition
          result <- run_evaluation("hOUwie", resolution, fixture,
                                   run_seed = run_seed)
          row <- data.frame(
            task_id = task_id,
            profile = profile_name,
            method = "hOUwie",
            tips = n_tips,
            data_seed = data_seed,
            run_seed = run_seed,
            evaluation_repetition = evaluation_repetition,
            resolution_name = "sampled histories",
            resolution = resolution,
            elapsed_seconds = as.numeric(result["elapsed_seconds"]),
            loglik = as.numeric(result["loglik"]),
            error = unname(result["error"]),
            stringsAsFactors = FALSE
          )
          append_result(row, evaluations_path)
        }
      }
    }

    for(resolution in profile$fitmultiou_eval_levels){
      for(evaluation_repetition in seq_len(profile$evaluation_repetitions)){
        task_id <- sprintf("%s-eval-fitmultiou-t%d-d%d-r%d-rep%d",
                           profile_name, n_tips, data_seed, resolution,
                           evaluation_repetition)
        if(!task_finished(task_id, evaluations_path)){
          cat("Running", task_id, "\n")
          run_seed <- data_seed * 10000L + evaluation_repetition
          result <- run_evaluation("fitmultiOU", resolution, fixture,
                                   run_seed = run_seed)
          row <- data.frame(
            task_id = task_id,
            profile = profile_name,
            method = "fitmultiOU",
            tips = n_tips,
            data_seed = data_seed,
            run_seed = run_seed,
            evaluation_repetition = evaluation_repetition,
            resolution_name = "discretization levels",
            resolution = resolution,
            elapsed_seconds = as.numeric(result["elapsed_seconds"]),
            loglik = as.numeric(result["loglik"]),
            error = unname(result["error"]),
            stringsAsFactors = FALSE
          )
          append_result(row, evaluations_path)
        }
      }
    }

    fit_methods <- c("hOUwie", "fitmultiOU")
    if(profile$include_parallel_houwie && cores > 1L){
      fit_methods <- c(fit_methods, "hOUwie parallel")
    }
    for(repetition in seq_len(profile$fit_repetitions)){
      for(method in fit_methods){
        method_id <- gsub("[^a-z0-9]+", "-", tolower(method))
        task_id <- sprintf("%s-fit-%s-t%d-d%d-rep%d", profile_name,
                           method_id, n_tips, data_seed, repetition)
        if(task_finished(task_id, fits_path)) next
        cat("Running", task_id, "\n")
        run_seed <- data_seed * 1000L + repetition
        result <- run_fit(method, fixture, run_seed)
        row <- data.frame(
          task_id = task_id,
          profile = profile_name,
          method = method,
          tips = n_tips,
          data_seed = data_seed,
          run_seed = run_seed,
          resolution_name = if(grepl("^hOUwie", method))
            "sampled histories" else "discretization levels",
          resolution = if(grepl("^hOUwie", method))
            profile$houwie_nsim else profile$fitmultiou_levels,
          starts = if(grepl("^hOUwie", method))
            profile$houwie_starts else 1L,
          cores = if(method == "hOUwie parallel")
            min(cores, profile$houwie_starts) else if(method == "hOUwie")
            1L else cores,
          max_evaluations = if(grepl("^hOUwie", method))
            profile$houwie_maxeval else profile$fitmultiou_maxit,
          elapsed_seconds = as.numeric(result["elapsed_seconds"]),
          loglik = as.numeric(result["loglik"]),
          theta_a = as.numeric(result["theta_a"]),
          theta_b = as.numeric(result["theta_b"]),
          alpha = as.numeric(result["alpha"]),
          sigma_sq = as.numeric(result["sigma_sq"]),
          q = as.numeric(result["q"]),
          evaluations = as.numeric(result["evaluations"]),
          convergence = as.numeric(result["convergence"]),
          theta_rmse = as.numeric(result["theta_rmse"]),
          alpha_log_error = as.numeric(result["alpha_log_error"]),
          sigma_log_error = as.numeric(result["sigma_log_error"]),
          q_log_error = as.numeric(result["q_log_error"]),
          error = unname(result["error"]),
          stringsAsFactors = FALSE
        )
        append_result(row, fits_path)
      }
    }
  }
}

cat("Finished profile", profile_name, "\n")
cat("Results:", results_dir, "\n")
