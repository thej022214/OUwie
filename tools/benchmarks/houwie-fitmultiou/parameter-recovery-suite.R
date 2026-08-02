#!/usr/bin/env Rscript

# A compact parameter-recovery sweep around the simulation in Liam Revell's
# 31 July 2026 post. Each hOUwie fit starts at the generating parameters, which
# makes this a deliberately favorable test of recovery rather than a search for
# a good initial point. Results are appended after every fit so the run can be
# resumed.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]),
                             mustWork = TRUE)
benchmark_dir <- dirname(script_file)
repo_root <- normalizePath(file.path(benchmark_dir, "..", "..", ".."),
                           mustWork = TRUE)

phytools_lib <- Sys.getenv("PHYTOOLS_LIB")
if(nzchar(phytools_lib)) .libPaths(c(phytools_lib, .libPaths()))

if(!requireNamespace("phytools", quietly = TRUE) ||
   !"fitmultiOU" %in% getNamespaceExports("phytools")){
  stop("A development phytools containing fitmultiOU is required.",
       call. = FALSE)
}
if(!requireNamespace("pkgload", quietly = TRUE)){
  stop("pkgload is required.", call. = FALSE)
}
pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                  helpers = FALSE)

scenarios <- data.frame(
  scenario = c(
    "baseline", "weak-alpha", "strong-alpha", "low-sigma-squared",
    "high-sigma-squared", "small-optimum-separation",
    "large-optimum-separation", "low-q", "high-q"
  ),
  theta_a = c(5, 5, 5, 5, 5, 6.5, 2.5, 5, 5),
  theta_b = c(10, 10, 10, 10, 10, 8.5, 12.5, 10, 10),
  alpha = c(0.5, 0.05, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  sigma_sq = c(0.2, 0.2, 0.2, 0.05, 0.8, 0.2, 0.2, 0.2, 0.2),
  q = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.05, 1),
  stringsAsFactors = FALSE
)

seed_values <- as.integer(strsplit(
  Sys.getenv("HOUWIE_SUITE_SEEDS", "1101,1102,1103"),
  ",", fixed = TRUE
)[[1]])
n_tips <- as.integer(Sys.getenv("HOUWIE_SUITE_TIPS", "50"))
n_sim <- as.integer(Sys.getenv("HOUWIE_SUITE_NSIM", "50"))
maxeval <- as.integer(Sys.getenv("HOUWIE_SUITE_MAXEVAL", "150"))
levels <- as.integer(Sys.getenv("HOUWIE_SUITE_LEVELS", "50"))
methods <- strsplit(
  Sys.getenv("HOUWIE_SUITE_METHODS", "hOUwie,fitmultiOU"),
  ",", fixed = TRUE
)[[1]]
if(any(!is.finite(c(seed_values, n_tips, n_sim, maxeval, levels))) ||
   any(c(n_tips, n_sim, maxeval, levels) < 1)){
  stop("Suite settings must be positive integers.", call. = FALSE)
}
if(!length(methods) || any(!methods %in% c("hOUwie", "fitmultiOU"))){
  stop("HOUWIE_SUITE_METHODS must contain hOUwie and/or fitmultiOU.",
       call. = FALSE)
}

results_dir <- file.path(benchmark_dir, "results", "parameter-recovery-suite")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(results_dir, "fits.csv")
if(identical(Sys.getenv("HOUWIE_BENCH_OVERWRITE"), "1")){
  unlink(results_path)
}

read_results <- function(){
  if(file.exists(results_path)) {
    read.csv(results_path, stringsAsFactors = FALSE)
  } else NULL
}

append_result <- function(row){
  old <- read_results()
  combined <- if(is.null(old)) row else rbind(old, row)
  write.csv(combined, results_path, row.names = FALSE)
}

simulate_fixture <- function(parameters, seed){
  states <- c("a", "b")
  for(attempt in 0:100){
    actual_seed <- seed + 100000L * attempt
    set.seed(actual_seed)
    phy <- phytools::pbtree(n = n_tips, scale = 10)
    Q <- matrix(parameters$q, 2, 2,
                dimnames = list(states, states))
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    root_state <- sample(states, 1)
    sim_tree <- phytools::sim.history(phy, Q, anc = root_state)
    discrete <- factor(phytools::getStates(sim_tree, "tips"),
                       levels = states)
    if(all(table(discrete) > 0)) break
  }
  if(any(table(discrete) == 0)){
    stop("Could not simulate both tip states for ", parameters$scenario,
         " seed ", seed, call. = FALSE)
  }
  alpha <- setNames(rep(parameters$alpha, 2), states)
  sigma_sq <- setNames(rep(parameters$sigma_sq, 2), states)
  theta <- setNames(c(parameters$theta_a, parameters$theta_b), states)
  trait <- phytools::multiOU(
    sim_tree, alpha, sigma_sq, theta, a0 = theta[root_state]
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
    root_state = root_state,
    actual_seed = actual_seed
  )
}

fit_houwie <- function(fixture, parameters, run_seed){
  set.seed(run_seed)
  hOUwie(
    phy = fixture$phy,
    data = fixture$data,
    rate.cat = 1,
    discrete_model = "ER",
    continuous_model = "OUM",
    nSim = n_sim,
    ip = unname(c(
      parameters$q, parameters$alpha, parameters$sigma_sq,
      parameters$theta_a, parameters$theta_b
    )),
    opts = list(
      algorithm = "NLOPT_LN_SBPLX",
      maxeval = maxeval,
      ftol_rel = .Machine$double.eps^0.25
    ),
    common_random_numbers = TRUE,
    n_starts = 1,
    ncores = 1,
    quiet = TRUE,
    sample_tips = FALSE,
    sample_nodes = FALSE,
    adaptive_sampling = FALSE
  )
}

fit_fitmultiou <- function(fixture, parameters){
  initial <- setNames(
    c(
      parameters$theta_a, parameters$theta_b, parameters$alpha,
      parameters$sigma_sq, parameters$q
    ),
    c("theta[a]", "theta[b]", "alpha", "sigsq", "q[1]")
  )
  phytools::fitmultiOU(
    fixture$phy,
    fixture$trait,
    fixture$discrete,
    model = "ER",
    levs = levels,
    parallel = FALSE,
    root = "mle",
    trace = 0,
    maxit = maxeval,
    rand_start = FALSE,
    init = initial
  )
}

extract_houwie <- function(fit){
  state_columns <- c("(1)", "(2)")
  c(
    loglik = fit$loglik,
    theta_a = fit$solution.cont["theta", state_columns[1]],
    theta_b = fit$solution.cont["theta", state_columns[2]],
    alpha = fit$solution.cont["alpha", state_columns[1]],
    sigma_sq = fit$solution.cont["sigma2", state_columns[1]],
    q = mean(fit$solution.disc[!is.na(fit$solution.disc)]),
    evaluations = if(is.null(fit$global_liks_mat)) NA_real_ else
      sum(fit$global_liks_mat[[1]] != 0)
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
    evaluations = unname(fit$opt_results$counts[["function"]])
  )
}

for(scenario_i in seq_len(nrow(scenarios))){
  parameters <- scenarios[scenario_i, ]
  for(data_seed in seed_values){
    fixture_seed <- data_seed + 1000L * scenario_i
    fixture <- simulate_fixture(parameters, fixture_seed)
    for(method in methods){
      task_id <- paste(parameters$scenario, data_seed, method, sep = "-")
      old <- read_results()
      if(!is.null(old) && task_id %in% old$task_id) next
      cat("Running", task_id, "\n")
      error_message <- NA_character_
      started <- proc.time()[["elapsed"]]
      fit <- tryCatch(
        if(method == "hOUwie"){
          fit_houwie(fixture, parameters, fixture_seed + 500000L)
        }else{
          fit_fitmultiou(fixture, parameters)
        },
        error = function(e){
          error_message <<- conditionMessage(e)
          NULL
        }
      )
      elapsed <- proc.time()[["elapsed"]] - started
      estimate <- if(is.null(fit)){
        setNames(rep(NA_real_, 7), c(
          "loglik", "theta_a", "theta_b", "alpha", "sigma_sq", "q",
          "evaluations"
        ))
      }else if(method == "hOUwie"){
        extract_houwie(fit)
      }else{
        extract_fitmultiou(fit)
      }
      row <- data.frame(
        task_id = task_id,
        scenario = parameters$scenario,
        data_seed = data_seed,
        actual_seed = fixture$actual_seed,
        method = method,
        tips = n_tips,
        resolution = if(method == "hOUwie") n_sim else levels,
        max_evaluations = maxeval,
        elapsed_seconds = elapsed,
        true_theta_a = parameters$theta_a,
        true_theta_b = parameters$theta_b,
        true_alpha = parameters$alpha,
        true_sigma_sq = parameters$sigma_sq,
        true_q = parameters$q,
        loglik = unname(estimate["loglik"]),
        theta_a = unname(estimate["theta_a"]),
        theta_b = unname(estimate["theta_b"]),
        alpha = unname(estimate["alpha"]),
        sigma_sq = unname(estimate["sigma_sq"]),
        q = unname(estimate["q"]),
        evaluations = unname(estimate["evaluations"]),
        error = error_message,
        stringsAsFactors = FALSE
      )
      append_result(row)
    }
  }
}

cat("Finished parameter recovery suite\n")
cat("Results:", results_path, "\n")
