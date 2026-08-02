#!/usr/bin/env Rscript

# Score the candidate parameter vectors from Liam Revell's 300-tip simulation
# with the deterministic mixture pruning, and compare the ranking against what
# hOUwie's sampled sum reports for the same vectors.

benchmark_dir <- normalizePath(dirname(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
))
repo_root <- normalizePath(file.path(benchmark_dir, "..", "..", ".."))
suppressPackageStartupMessages({
  library(ape)
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

candidates <- data.frame(
  candidate = c("Generating parameters", "fitmultiOU estimate from post",
                "hOUwie estimate from post"),
  theta_a = c(5, 5.0209, 4.4566072),
  theta_b = c(10, 10.1985, 88.5531307),
  alpha = c(0.5, 0.4919, 0.0123027),
  sigma_sq = c(0.2, 0.1839, 0.7764039),
  q = c(0.3, 0.347199, 0.1187636),
  stringsAsFactors = FALSE
)

caps <- c(2, 4, 8)
resolutions <- c(1L, 2L, 4L)
rows <- list()
for(i in seq_len(nrow(candidates))){
  candidate <- candidates[i, ]
  Q <- matrix(c(-candidate$q, candidate$q, candidate$q, -candidate$q), 2, 2,
              byrow = TRUE)
  alpha <- rep(candidate$alpha, 2)
  sigma_sq <- rep(candidate$sigma_sq, 2)
  theta <- c(candidate$theta_a, candidate$theta_b)
  for(resolution in resolutions){
    for(cap in caps){
      started <- Sys.time()
      value <- mixture_pruning_loglik(phy, tip_states, tip_values, Q, alpha,
                                      sigma_sq, theta, root_p = c(0.5, 0.5),
                                      resolution = resolution,
                                      max_components = cap)
      rows[[length(rows) + 1L]] <- data.frame(
        candidate = candidate$candidate, resolution = resolution,
        max_components = cap, loglik = value,
        seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
      cat(sprintf("%-30s res=%d cap=%2d  lnL=%12.4f  %5.1fs\n",
                  candidate$candidate, resolution, cap, value,
                  rows[[length(rows)]]$seconds))
    }
  }
}
result <- do.call(rbind, rows)
dir.create(file.path(benchmark_dir, "results", "mixture-pruning"),
           recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(benchmark_dir, "results", "mixture-pruning",
                            "candidate-likelihoods.csv"), row.names = FALSE)
print(reshape(result[, c("candidate", "resolution", "max_components",
                         "loglik")],
              idvar = c("resolution", "max_components"),
              timevar = "candidate", direction = "wide"),
      row.names = FALSE, digits = 8)
