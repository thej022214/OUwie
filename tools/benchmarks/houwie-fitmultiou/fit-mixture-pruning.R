#!/usr/bin/env Rscript

# Optimize the deterministic mixture-pruning likelihood on the 300-tip fixture
# and report where the MLE lands relative to the generating parameters. This is
# the test the likelihood-only comparison cannot make: scoring three candidate
# vectors correctly does not prove that searching the surface recovers them.

benchmark_dir <- normalizePath(dirname(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
))
repo_root <- normalizePath(file.path(benchmark_dir, "..", "..", ".."))
suppressPackageStartupMessages({
  library(ape)
  library(nloptr)
  pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                    helpers = FALSE)
})
source(file.path(benchmark_dir, "mixture-pruning-prototype.R"))

fixture <- readRDS(file.path(benchmark_dir, "results", "scaling",
                             "fixture-tips300-seed8.rds"))
phy <- ape::reorder.phylo(fixture$phy, "pruningwise")
tip_states <- as.integer(fixture$discrete)
names(tip_states) <- names(fixture$discrete)
tip_values <- fixture$trait
Tmax <- max(ape::branching.times(phy))

# parameters are c(q, alpha, sigma_sq, theta_a, theta_b), searched on the log
# scale inside hOUwie's own bounds
truth <- c(q = 0.3, alpha = 0.5, sigma_sq = 0.2, theta_a = 5, theta_b = 10)
lower <- log(c(1 / (Tmax * 10000), 1e-10, 1e-10,
               rep(min(tip_values) / 10, 2)))
upper <- log(c(1 / (Tmax * 0.0001), log(2) / (0.01 * Tmax),
               log(2) / (0.01 * Tmax), rep(max(tip_values) * 10, 2)))

objective_for <- function(resolution, max_components){
  function(p){
    value <- exp(p)
    Q <- matrix(c(-value[1], value[1], value[1], -value[1]), 2, 2,
                byrow = TRUE)
    result <- try(mixture_pruning_loglik(
      phy, tip_states, tip_values, Q,
      alpha = rep(value[2], 2), sigma_sq = rep(value[3], 2),
      theta = c(value[4], value[5]), root_p = c(0.5, 0.5),
      resolution = resolution, max_components = max_components
    ), silent = TRUE)
    if(inherits(result, "try-error") || !is.finite(result)) return(1e10)
    -result
  }
}

# hOUwie's own default start, plus perturbations of it, so the search is not
# handed the answer
set.seed(4021)
base_start <- c(10 / sum(phy$edge.length), log(2) / Tmax, log(2) / Tmax,
                unname(quantile(tip_values, 0.25)),
                unname(quantile(tip_values, 0.75)))
starts <- c(list(base_start), lapply(1:3, function(i){
  base_start * exp(stats::rnorm(5, 0, 0.4))
}))

run_stage <- function(resolution, max_components, start_list, label){
  rows <- list()
  for(i in seq_along(start_list)){
    start <- pmin(pmax(log(start_list[[i]]), lower + 1e-6), upper - 1e-6)
    began <- Sys.time()
    fit <- nloptr(x0 = start, eval_f = objective_for(resolution,
                                                     max_components),
                  lb = lower, ub = upper,
                  opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1500,
                              ftol_rel = .Machine$double.eps^0.25))
    elapsed <- as.numeric(difftime(Sys.time(), began, units = "mins"))
    estimate <- exp(fit$solution)
    rows[[i]] <- data.frame(
      stage = label, resolution = resolution,
      max_components = max_components, start = i,
      loglik = -fit$objective, q = estimate[1], alpha = estimate[2],
      sigma_sq = estimate[3], theta_a = estimate[4], theta_b = estimate[5],
      evaluations = fit$iterations, minutes = elapsed
    )
    cat(sprintf(
      "%-8s start %d  lnL=%9.3f  q=%.4f a=%.4f s2=%.4f th=(%.3f, %.3f)  %.1f min\n",
      label, i, -fit$objective, estimate[1], estimate[2], estimate[3],
      estimate[4], estimate[5], elapsed))
    utils::flush.console()
  }
  do.call(rbind, rows)
}

# coarse settings first, because they are the ones fast enough to search from
# several starts; then refine the winner where the approximation is tighter
coarse <- run_stage(2L, 2L, starts, "coarse")
best <- coarse[which.max(coarse$loglik), ]
refined <- run_stage(4L, 4L, list(unlist(best[, c("q", "alpha", "sigma_sq",
                                                  "theta_a", "theta_b")])),
                     "refined")

result <- rbind(coarse, refined)
out_dir <- file.path(benchmark_dir, "results", "mixture-pruning")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(out_dir, "pruning-mle.csv"), row.names = FALSE)

cat("\ngenerating parameters:\n")
print(truth)
cat("\nall fits:\n")
print(result[, c("stage", "start", "loglik", "q", "alpha", "sigma_sq",
                 "theta_a", "theta_b")], row.names = FALSE, digits = 4)
