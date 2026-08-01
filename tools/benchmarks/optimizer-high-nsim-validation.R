#!/usr/bin/env Rscript

# Re-score the completed optimizer comparison at a larger, fixed nSim. This is
# validation only: it does not alter or resume any optimization.

benchmark_library <- Sys.getenv(
  "HOUWIE_BENCH_LIB",
  unset = "/private/tmp/houwie-benchmark-lib"
)
if (dir.exists(benchmark_library)) {
  .libPaths(c(benchmark_library, .libPaths()))
}

devtools::load_all(quiet = TRUE)
data.table::setDTthreads(1L)

results_dir <- Sys.getenv(
  "HOUWIE_BENCH_RESULTS",
  unset = file.path(
    normalizePath("tools/benchmarks", mustWork = TRUE),
    "optimizer-comparison-results-moderate"
  )
)
input <- readRDS(file.path(results_dir, "benchmark-input.rds"))
fits <- readRDS(file.path(results_dir, "validated-fits.rds"))
high_nsim <- suppressWarnings(as.integer(Sys.getenv("HOUWIE_HIGH_NSIM", "250")))
n_validation <- suppressWarnings(as.integer(Sys.getenv("HOUWIE_HIGH_NVALIDATION", "10")))
cores <- suppressWarnings(as.integer(Sys.getenv("HOUWIE_BENCH_CORES", "3")))
if (!is.finite(high_nsim) || high_nsim < 1L) high_nsim <- 250L
if (!is.finite(n_validation) || n_validation < 1L) n_validation <- 10L
if (!is.finite(cores) || cores < 1L) cores <- 1L

validation_seeds <- input$validation_seeds[seq_len(n_validation)]
simulated <- input$simulated

candidates <- lapply(fits, function(fit) {
  list(
    optimizer = fit$optimizer,
    start_id = fit$start_id,
    parameters = exp(fit$solution)
  )
})
candidates[[length(candidates) + 1L]] <- list(
  optimizer = "TRUTH",
  start_id = 0L,
  parameters = unname(simulated$truth)
)

checkpoint_dir <- file.path(results_dir, paste0("validation-nsim-", high_nsim))
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

score_candidate <- function(candidate) {
  checkpoint <- file.path(
    checkpoint_dir,
    paste0(tolower(candidate$optimizer), "-start-", candidate$start_id, ".rds")
  )
  if (file.exists(checkpoint)) return(readRDS(checkpoint))

  loglik <- vapply(validation_seeds, function(seed) {
    set.seed(seed)
    fit <- tryCatch(
      suppressWarnings(hOUwie(
        phy = simulated$phy,
        data = simulated$data,
        rate.cat = 2L,
        discrete_model = "ER",
        continuous_model = "OUM",
        null.model = TRUE,
        nSim = high_nsim,
        p = candidate$parameters,
        quiet = TRUE,
        sample_tips = FALSE,
        sample_nodes = FALSE,
        adaptive_sampling = FALSE
      )),
      error = function(error) NULL
    )
    if (is.null(fit)) NA_real_ else fit$loglik
  }, numeric(1))

  result <- c(candidate, list(loglik = loglik))
  saveRDS(result, checkpoint)
  result
}

scores <- parallel::mclapply(
  candidates,
  score_candidate,
  mc.cores = min(cores, length(candidates)),
  mc.preschedule = FALSE
)

score_rows <- do.call(rbind, lapply(scores, function(score) {
  data.frame(
    optimizer = score$optimizer,
    start_id = score$start_id,
    validation_id = seq_along(score$loglik),
    loglik = score$loglik,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  score_rows,
  file.path(results_dir, paste0("validation-nsim-", high_nsim, ".csv")),
  row.names = FALSE
)

candidate_summary <- do.call(rbind, lapply(scores, function(score) {
  data.frame(
    optimizer = score$optimizer,
    start_id = score$start_id,
    mean_loglik = mean(score$loglik, na.rm = TRUE),
    median_loglik = stats::median(score$loglik, na.rm = TRUE),
    sd_loglik = stats::sd(score$loglik, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  candidate_summary,
  file.path(results_dir, paste0("validation-nsim-", high_nsim, "-summary.csv")),
  row.names = FALSE
)

optimizer_summary <- do.call(
  rbind,
  lapply(split(candidate_summary, candidate_summary$optimizer), function(rows) {
    data.frame(
      optimizer = rows$optimizer[1L],
      mean_loglik = mean(rows$mean_loglik),
      median_loglik = stats::median(rows$mean_loglik),
      between_start_sd = if (nrow(rows) > 1L) stats::sd(rows$mean_loglik) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
)
print(optimizer_summary, row.names = FALSE)
