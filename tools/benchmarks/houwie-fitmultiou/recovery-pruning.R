#!/usr/bin/env Rscript

# Parameter recovery for the deterministic mixture pruning, on exactly the
# fixtures parameter-recovery-suite.R uses. Two arms are run:
#
#   pruning  - the mixture pruning likelihood
#   hOUwie   - the sampled likelihood, re-run after the map.states fix
#
# so the table separates what the labelling bug cost from what the sampled
# estimator costs. Both start at the generating parameters, which makes this a
# test of whether the surface holds the truth rather than of whether a search
# can find it. Results are appended after each fit so the run can be resumed.

script_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
benchmark_dir <- dirname(script_file)
repo_root <- normalizePath(file.path(benchmark_dir, "..", "..", ".."))
suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(nloptr)
  pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                    helpers = FALSE)
})
source(file.path(benchmark_dir, "mixture-pruning-prototype.R"))

scenarios <- data.frame(
  scenario = c("baseline", "weak-alpha", "strong-alpha", "low-sigma-squared",
               "high-sigma-squared", "small-optimum-separation",
               "large-optimum-separation", "low-q", "high-q"),
  theta_a = c(5, 5, 5, 5, 5, 6.5, 2.5, 5, 5),
  theta_b = c(10, 10, 10, 10, 10, 8.5, 12.5, 10, 10),
  alpha = c(0.5, 0.05, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  sigma_sq = c(0.2, 0.2, 0.2, 0.05, 0.8, 0.2, 0.2, 0.2, 0.2),
  q = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.05, 1),
  stringsAsFactors = FALSE
)
seed_values <- c(1101L, 1102L, 1103L)
n_tips <- 50L
n_sim <- 50L
maxeval <- 150L
resolution <- 4L
max_components <- 4L

results_dir <- file.path(benchmark_dir, "results", "recovery-pruning")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(results_dir, "fits.csv")
read_results <- function(){
  if(file.exists(results_path)) read.csv(results_path, stringsAsFactors = FALSE)
  else NULL
}

# identical to parameter-recovery-suite.R, so the fixtures match case for case
simulate_fixture <- function(parameters, seed){
  states <- c("a", "b")
  for(attempt in 0:100){
    actual_seed <- seed + 100000L * attempt
    set.seed(actual_seed)
    phy <- phytools::pbtree(n = n_tips, scale = 10)
    Q <- matrix(parameters$q, 2, 2, dimnames = list(states, states))
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    root_state <- sample(states, 1)
    sim_tree <- phytools::sim.history(phy, Q, anc = root_state)
    discrete <- factor(phytools::getStates(sim_tree, "tips"), levels = states)
    if(all(table(discrete) > 0)) break
  }
  alpha <- setNames(rep(parameters$alpha, 2), states)
  sigma_sq <- setNames(rep(parameters$sigma_sq, 2), states)
  theta <- setNames(c(parameters$theta_a, parameters$theta_b), states)
  trait <- phytools::multiOU(sim_tree, alpha, sigma_sq, theta,
                             a0 = theta[root_state])
  list(phy = phy, trait = trait, discrete = discrete,
       data = data.frame(species = names(trait), state = discrete,
                         trait = as.numeric(trait), stringsAsFactors = FALSE),
       actual_seed = actual_seed)
}

fit_pruning <- function(fixture, parameters){
  phy <- ape::reorder.phylo(fixture$phy, "pruningwise")
  tip_states <- as.integer(fixture$discrete)
  names(tip_states) <- names(fixture$trait)
  Tmax <- max(ape::branching.times(phy))
  lower <- log(c(1 / (Tmax * 10000), 1e-10, 1e-10,
                 rep(min(fixture$trait) / 10, 2)))
  upper <- log(c(1 / (Tmax * 0.0001), log(2) / (0.01 * Tmax),
                 log(2) / (0.01 * Tmax), rep(max(fixture$trait) * 10, 2)))
  objective <- function(p){
    v <- exp(p)
    Q <- matrix(c(-v[1], v[1], v[1], -v[1]), 2, 2, byrow = TRUE)
    out <- try(mixture_pruning_loglik(
      phy, tip_states, fixture$trait, Q, alpha = rep(v[2], 2),
      sigma_sq = rep(v[3], 2), theta = c(v[4], v[5]), root_p = c(0.5, 0.5),
      resolution = resolution, max_components = max_components), silent = TRUE)
    if(inherits(out, "try-error") || !is.finite(out)) return(1e10)
    -as.numeric(out)
  }
  start <- log(c(parameters$q, parameters$alpha, parameters$sigma_sq,
                 parameters$theta_a, parameters$theta_b))
  start <- pmin(pmax(start, lower + 1e-6), upper - 1e-6)
  fit <- nloptr(x0 = start, eval_f = objective, lb = lower, ub = upper,
                opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = maxeval,
                            ftol_rel = .Machine$double.eps^0.25))
  estimate <- exp(fit$solution)
  c(loglik = -fit$objective, theta_a = estimate[4], theta_b = estimate[5],
    alpha = estimate[2], sigma_sq = estimate[3], q = estimate[1])
}

fit_houwie <- function(fixture, parameters, run_seed){
  set.seed(run_seed)
  fit <- hOUwie(phy = fixture$phy, data = fixture$data, rate.cat = 1,
                discrete_model = "ER", continuous_model = "OUM", nSim = n_sim,
                ip = unname(c(parameters$q, parameters$alpha,
                              parameters$sigma_sq, parameters$theta_a,
                              parameters$theta_b)),
                opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = maxeval,
                            ftol_rel = .Machine$double.eps^0.25),
                common_random_numbers = TRUE, n_starts = 1, ncores = 1,
                quiet = TRUE, sample_tips = FALSE, sample_nodes = FALSE,
                adaptive_sampling = FALSE)
  c(loglik = fit$loglik,
    theta_a = fit$solution.cont["theta", "(1)"],
    theta_b = fit$solution.cont["theta", "(2)"],
    alpha = fit$solution.cont["alpha", "(1)"],
    sigma_sq = fit$solution.cont["sigma2", "(1)"],
    q = mean(fit$solution.disc[!is.na(fit$solution.disc)]))
}

for(scenario_i in seq_len(nrow(scenarios))){
  parameters <- scenarios[scenario_i, ]
  for(data_seed in seed_values){
    fixture <- simulate_fixture(parameters, data_seed + 1000L * scenario_i)
    for(method in c("pruning", "hOUwie")){
      task_id <- paste(parameters$scenario, data_seed, method, sep = "-")
      old <- read_results()
      if(!is.null(old) && task_id %in% old$task_id) next
      error_message <- NA_character_
      started <- proc.time()[["elapsed"]]
      estimate <- tryCatch(
        if(method == "pruning") fit_pruning(fixture, parameters)
        else fit_houwie(fixture, parameters,
                        data_seed + 1000L * scenario_i + 500000L),
        error = function(e){ error_message <<- conditionMessage(e); NULL })
      if(is.null(estimate)){
        estimate <- setNames(rep(NA_real_, 6),
          c("loglik", "theta_a", "theta_b", "alpha", "sigma_sq", "q"))
      }
      row <- data.frame(
        task_id = task_id, scenario = parameters$scenario,
        data_seed = data_seed, method = method,
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        true_theta_a = parameters$theta_a, true_theta_b = parameters$theta_b,
        true_alpha = parameters$alpha, true_sigma_sq = parameters$sigma_sq,
        true_q = parameters$q,
        loglik = unname(estimate["loglik"]),
        theta_a = unname(estimate["theta_a"]),
        theta_b = unname(estimate["theta_b"]),
        alpha = unname(estimate["alpha"]),
        sigma_sq = unname(estimate["sigma_sq"]),
        q = unname(estimate["q"]),
        error = error_message, stringsAsFactors = FALSE)
      write.csv(if(is.null(old)) row else rbind(old, row), results_path,
                row.names = FALSE)
      cat(sprintf("%-34s a=%7.3f s2=%7.3f q=%6.3f th=(%6.2f,%6.2f) %5.1fs\n",
                  task_id, row$alpha, row$sigma_sq, row$q, row$theta_a,
                  row$theta_b, row$elapsed_seconds))
      utils::flush.console()
    }
  }
}
cat("done ->", results_path, "\n")
