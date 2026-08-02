#!/usr/bin/env Rscript

# Evaluate three biologically meaningful parameter vectors on the exact
# 300-tip simulation from Liam Revell's 31 July 2026 post. The script compares
# each method's likelihood ranking, repeats the hOUwie calculation on shared
# history banks, and checks its continuous likelihood on the known generating
# history. Absolute likelihood scales are not compared across methods.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
benchmark_dir <- dirname(script_file)
repo_root <- normalizePath(file.path(benchmark_dir, "..", "..", ".."),
                           mustWork = TRUE)

phytools_lib <- Sys.getenv("PHYTOOLS_LIB")
if(nzchar(phytools_lib)){
  if(!dir.exists(phytools_lib)){
    stop("PHYTOOLS_LIB does not exist: ", phytools_lib, call. = FALSE)
  }
  .libPaths(c(normalizePath(phytools_lib), .libPaths()))
}
if(!requireNamespace("phytools", quietly = TRUE) ||
   !"fitmultiOU" %in% getNamespaceExports("phytools")){
  stop(
    "The development version of phytools containing fitmultiOU is required.",
    call. = FALSE
  )
}
if(!requireNamespace("pkgload", quietly = TRUE)){
  stop("pkgload is required to load the current local OUwie source.",
       call. = FALSE)
}
pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                  helpers = FALSE)

fixture_path <- file.path(
  benchmark_dir, "results", "scaling", "fixture-tips300-seed8.rds"
)
if(!file.exists(fixture_path)){
  stop(
    "The exact 300-tip fixture is missing. Run the scaling profile first.",
    call. = FALSE
  )
}
fixture <- readRDS(fixture_path)

candidates <- data.frame(
  candidate = c(
    "Generating parameters",
    "fitmultiOU estimate from post",
    "hOUwie estimate from post"
  ),
  theta_a = c(5, 5.0209, 4.4566072),
  theta_b = c(10, 10.1985, 88.5531307),
  alpha = c(0.5, 0.4919, 0.0123027),
  sigma_sq = c(0.2, 0.1839, 0.7764039),
  q = c(0.3, 0.347199, 0.1187636),
  stringsAsFactors = FALSE
)

repetitions <- suppressWarnings(as.integer(
  Sys.getenv("HOUWIE_CANDIDATE_REPS", "10")
))
n_sim <- suppressWarnings(as.integer(
  Sys.getenv("HOUWIE_CANDIDATE_NSIM", "100")
))
if(is.na(repetitions) || repetitions < 1L || is.na(n_sim) || n_sim < 1L){
  stop("HOUWIE_CANDIDATE_REPS and HOUWIE_CANDIDATE_NSIM must be positive integers.",
       call. = FALSE)
}

results_dir <- file.path(benchmark_dir, "results", "parameter-comparison")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(results_dir, "candidate-likelihoods.csv")
candidates_path <- file.path(results_dir, "candidate-parameters.csv")
fixed_results_path <- file.path(
  results_dir, "fixed-history-candidate-likelihoods.csv"
)
known_history_path <- file.path(
  results_dir, "known-history-continuous-likelihoods.csv"
)
write.csv(candidates, candidates_path, row.names = FALSE)

overwrite <- identical(Sys.getenv("HOUWIE_BENCH_OVERWRITE"), "1")
if(overwrite){
  unlink(c(results_path, fixed_results_path, known_history_path))
}

read_results <- function(){
  if(file.exists(results_path))
    read.csv(results_path, stringsAsFactors = FALSE) else NULL
}

task_finished <- function(task_id){
  old <- read_results()
  !is.null(old) && task_id %in% old$task_id
}

append_result <- function(row){
  old <- read_results()
  combined <- if(is.null(old)) row else rbind(old, row)
  write.csv(combined, results_path, row.names = FALSE)
}

for(repetition in seq_len(repetitions)){
  run_seed <- 90000L + repetition
  for(i in seq_len(nrow(candidates))){
    candidate <- candidates[i, ]
    candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidate$candidate))
    task_id <- sprintf("houwie-%s-rep%d", candidate_id, repetition)
    if(task_finished(task_id)) next
    cat("Running", task_id, "\n")
    set.seed(run_seed)
    started <- proc.time()[["elapsed"]]
    error_message <- NA_character_
    fit <- tryCatch(
      hOUwie(
        phy = fixture$phy,
        data = fixture$data,
        rate.cat = 1,
        discrete_model = "ER",
        continuous_model = "OUM",
        nSim = n_sim,
        p = unname(c(
          candidate$q,
          candidate$alpha,
          candidate$sigma_sq,
          candidate$theta_a,
          candidate$theta_b
        )),
        common_random_numbers = TRUE,
        quiet = TRUE,
        sample_tips = FALSE,
        sample_nodes = FALSE,
        adaptive_sampling = FALSE
      ),
      error = function(e){
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    elapsed <- proc.time()[["elapsed"]] - started
    row <- data.frame(
      task_id = task_id,
      method = "hOUwie",
      candidate = candidate$candidate,
      repetition = repetition,
      run_seed = run_seed,
      resolution = n_sim,
      elapsed_seconds = elapsed,
      loglik = if(is.null(fit)) NA_real_ else fit$loglik,
      discrete_loglik = if(is.null(fit)) NA_real_ else fit$DiscLik,
      continuous_loglik = if(is.null(fit)) NA_real_ else fit$ContLik,
      error = error_message,
      stringsAsFactors = FALSE
    )
    append_result(row)
  }
}

# fitmultiOU is deterministic for a fixed grid and parameter vector, so a single
# evaluation of each candidate is sufficient. Building the closure is included
# in the first task only; elapsed time is retained solely for diagnostics.
fitmulti_setup <- NULL
for(i in seq_len(nrow(candidates))){
  candidate <- candidates[i, ]
  candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidate$candidate))
  task_id <- sprintf("fitmultiou-%s", candidate_id)
  if(task_finished(task_id)) next
  cat("Running", task_id, "\n")
  started <- proc.time()[["elapsed"]]
  error_message <- NA_character_
  loglik <- tryCatch(
    {
      if(is.null(fitmulti_setup)){
        initial <- setNames(
          c(5, 10, 0.5, 0.2, 0.3),
          c("theta[a]", "theta[b]", "alpha", "sigsq", "q[1]")
        )
        fitmulti_setup <- phytools::fitmultiOU(
          fixture$phy,
          fixture$trait,
          fixture$discrete,
          model = "ER",
          levs = 100,
          parallel = FALSE,
          root = "mle",
          trace = 0,
          maxit = 0,
          rand_start = FALSE,
          init = initial
        )
      }
      fitmulti_setup$lik(
        theta = unname(c(candidate$theta_a, candidate$theta_b)),
        alpha = log(candidate$alpha),
        sigsq = log(candidate$sigma_sq),
        q = log(candidate$q)
      )
    },
    error = function(e){
      error_message <<- conditionMessage(e)
      NA_real_
    }
  )
  elapsed <- proc.time()[["elapsed"]] - started
  row <- data.frame(
    task_id = task_id,
    method = "fitmultiOU",
    candidate = candidate$candidate,
    repetition = 1L,
    run_seed = NA_integer_,
    resolution = 100L,
    elapsed_seconds = elapsed,
    loglik = loglik,
    discrete_loglik = NA_real_,
    continuous_loglik = NA_real_,
    error = error_message,
    stringsAsFactors = FALSE
  )
  append_result(row)
}

# Separate the OU calculation from the uncertainty in the discrete history by
# evaluating every parameter vector on the exact history that generated the
# traits. This history is available only because the dataset is simulated.
if(!file.exists(known_history_path)){
  observed_traits <- setNames(c("a", "b"), c("1", "2"))
  known_history <- OUwie:::makeMapEdgesNumeric(
    fixture$sim_tree, observed_traits
  )
  data_ou <- OUwie:::organizeHOUwieDat(
    fixture$data, "none", TRUE
  )$data.ou
  all_paths <- lapply(
    seq_len(ape::Ntip(known_history) + ape::Nnode(known_history)),
    function(node) OUwie:::getPathToRoot(known_history, node)
  )
  tip_paths <- all_paths[seq_len(ape::Ntip(known_history))]
  known_rows <- lapply(seq_len(nrow(candidates)), function(i){
    candidate <- candidates[i, ]
    continuous_loglik <- OUwie:::OUwie.basic(
      known_history,
      data_ou,
      simmap.tree = TRUE,
      scaleHeight = FALSE,
      alpha = rep(candidate$alpha, 2),
      sigma.sq = rep(candidate$sigma_sq, 2),
      theta = unname(c(candidate$theta_a, candidate$theta_b)),
      algorithm = "three.point",
      tip.paths = tip_paths,
      tip.fog = "none"
    )
    data.frame(
      candidate = candidate$candidate,
      continuous_loglik = continuous_loglik,
      stringsAsFactors = FALSE
    )
  })
  write.csv(do.call(rbind, known_rows), known_history_path, row.names = FALSE)
}

# hOUwie normally samples a fresh history bank at every parameter vector. To
# hold that source of variation still, sample one bank under each candidate and
# evaluate all three candidates on every bank. The maps returned by hOUwie are
# already numbered 1 and 2, so hOUwie.fixed must not translate their labels a
# second time.
if(!file.exists(fixed_results_path)){
  fixed_rows <- list()
  row_index <- 0L
  for(bank_i in seq_len(nrow(candidates))){
    bank_source <- candidates[bank_i, ]
    cat("Sampling fixed history bank under", bank_source$candidate, "\n")
    set.seed(92000L + bank_i)
    bank_fit <- hOUwie(
      phy = fixture$phy,
      data = fixture$data,
      rate.cat = 1,
      discrete_model = "ER",
      continuous_model = "OUM",
      nSim = n_sim,
      p = unname(c(
        bank_source$q,
        bank_source$alpha,
        bank_source$sigma_sq,
        bank_source$theta_a,
        bank_source$theta_b
      )),
      common_random_numbers = TRUE,
      quiet = TRUE,
      sample_tips = FALSE,
      sample_nodes = FALSE,
      adaptive_sampling = FALSE
    )
    for(candidate_i in seq_len(nrow(candidates))){
      candidate <- candidates[candidate_i, ]
      cat("Evaluating", candidate$candidate, "on that fixed bank\n")
      error_message <- NA_character_
      fixed_fit <- tryCatch(
        hOUwie.fixed(
          simmaps = bank_fit$simmaps,
          data = fixture$data,
          rate.cat = 1,
          discrete_model = "ER",
          continuous_model = "OUM",
          p = unname(c(
            candidate$q,
            candidate$alpha,
            candidate$sigma_sq,
            candidate$theta_a,
            candidate$theta_b
          )),
          make_numeric = FALSE,
          quiet = TRUE,
          sample_tips = FALSE,
          sample_nodes = FALSE,
          adaptive_sampling = FALSE
        ),
        error = function(e){
          error_message <<- conditionMessage(e)
          NULL
        }
      )
      row_index <- row_index + 1L
      fixed_rows[[row_index]] <- data.frame(
        bank_source = bank_source$candidate,
        candidate = candidate$candidate,
        resolution = n_sim,
        loglik = if(is.null(fixed_fit)) NA_real_ else fixed_fit$loglik,
        discrete_loglik = if(is.null(fixed_fit)) NA_real_ else fixed_fit$DiscLik,
        continuous_loglik = if(is.null(fixed_fit)) NA_real_ else fixed_fit$ContLik,
        error = error_message,
        stringsAsFactors = FALSE
      )
    }
  }
  write.csv(do.call(rbind, fixed_rows), fixed_results_path, row.names = FALSE)
}

cat("Finished parameter comparison\n")
cat("Results:", results_path, "\n")
