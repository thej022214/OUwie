#!/usr/bin/env Rscript

# Score a recovery table against the four failure signatures reported for the
# original 27 hOUwie runs: alpha below truth, q below truth, sigma squared above
# truth, and the two optima estimated in the wrong order.

script_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
path <- file.path(dirname(script_file), "results", "recovery-pruning",
                  "fits.csv")
fits <- read.csv(path, stringsAsFactors = FALSE)
fits <- fits[is.na(fits$error), ]

signatures <- function(block){
  ordering_true <- sign(block$true_theta_b - block$true_theta_a)
  data.frame(
    method = block$method[1],
    n = nrow(block),
    alpha_below = sum(block$alpha < block$true_alpha),
    q_below = sum(block$q < block$true_q),
    sigma_above = sum(block$sigma_sq > block$true_sigma_sq),
    optima_reversed = sum(sign(block$theta_b - block$theta_a) != ordering_true),
    median_alpha_ratio = median(block$alpha / block$true_alpha),
    median_q_ratio = median(block$q / block$true_q),
    median_sigma_ratio = median(block$sigma_sq / block$true_sigma_sq),
    median_theta_gap_ratio = median(
      abs(block$theta_b - block$theta_a) /
        abs(block$true_theta_b - block$true_theta_a)),
    median_seconds = median(block$elapsed_seconds),
    stringsAsFactors = FALSE
  )
}

cat("failure signatures, counts out of n (lower is better for the first four)\n")
print(do.call(rbind, lapply(split(fits, fits$method), signatures)),
      row.names = FALSE, digits = 3)

cat("\nmedian |estimate / truth - 1| by scenario and method\n")
relative <- function(block){
  data.frame(scenario = block$scenario[1], method = block$method[1],
             alpha = median(abs(block$alpha / block$true_alpha - 1)),
             sigma_sq = median(abs(block$sigma_sq / block$true_sigma_sq - 1)),
             q = median(abs(block$q / block$true_q - 1)),
             theta_gap = median(abs(
               (block$theta_b - block$theta_a) /
                 (block$true_theta_b - block$true_theta_a) - 1)),
             stringsAsFactors = FALSE)
}
summary_table <- do.call(rbind, lapply(
  split(fits, list(fits$scenario, fits$method), drop = TRUE), relative))
print(summary_table[order(summary_table$scenario, summary_table$method), ],
      row.names = FALSE, digits = 2)
