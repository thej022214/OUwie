#!/usr/bin/env Rscript

# Compare derivative-free optimizers on a deliberately difficult simulated
# hOUwie problem. The latent maps remain parameter-dependent. During each fit,
# common random numbers make the Monte Carlo objective reproducible; all fits
# are then scored on the same fresh, independent random-number streams.

benchmark_library <- Sys.getenv(
  "HOUWIE_BENCH_LIB",
  unset = "/private/tmp/houwie-benchmark-lib"
)
if (dir.exists(benchmark_library)) {
  .libPaths(c(benchmark_library, .libPaths()))
}

required_packages <- c("devtools", "nloptr", "phytools", "crs")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing benchmark packages: ", paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

devtools::load_all(quiet = TRUE)
data.table::setDTthreads(1L)

benchmark_dir <- normalizePath("tools/benchmarks", mustWork = TRUE)
results_dir <- Sys.getenv(
  "HOUWIE_BENCH_RESULTS",
  unset = file.path(benchmark_dir, "optimizer-comparison-results-moderate")
)
fits_dir <- file.path(results_dir, "fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

integer_setting <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!is.finite(value) || value < 1L) default else value
}

n_sim <- integer_setting("HOUWIE_BENCH_NSIM", 25L)
n_starts <- integer_setting("HOUWIE_BENCH_NSTARTS", 4L)
n_validation <- integer_setting("HOUWIE_BENCH_NVALIDATION", 30L)
n_tips <- integer_setting("HOUWIE_BENCH_NTIPS", 96L)
tolerance <- .Machine$double.eps^0.25
benchmark_seed <- 20260731L
optimizer_names <- strsplit(
  Sys.getenv("HOUWIE_BENCH_OPTIMIZERS", "SBPLX,COBYLA,NOMAD"),
  ",",
  fixed = TRUE
)[[1L]]
optimizer_names <- trimws(optimizer_names)

simulate_challenging_data <- function(seed = benchmark_seed, n_tips = 96L) {
  set.seed(seed)
  phy <- phytools::pbtree(n = n_tips, scale = 1)

  state_names <- c("1A", "2A", "1B", "2B")
  q_matrix <- matrix(
    0,
    nrow = 4L,
    ncol = 4L,
    dimnames = list(state_names, state_names)
  )
  q_matrix["1A", "2A"] <- q_matrix["2A", "1A"] <- 0.8
  q_matrix["1B", "2B"] <- q_matrix["2B", "1B"] <- 1.5
  q_matrix["1A", "1B"] <- q_matrix["2A", "2B"] <- 0.3
  q_matrix["1B", "1A"] <- q_matrix["2B", "2A"] <- 0.2
  diag(q_matrix) <- -rowSums(q_matrix)

  set.seed(seed + 14L)
  full_history <- phytools::sim.history(
    phy,
    q_matrix,
    anc = NULL,
    direction = "row_to_column"
  )

  full_tip_data <- data.frame(
    species = phy$tip.label,
    state = unname(full_history$states[phy$tip.label])
  )
  set.seed(seed + 1L)
  simulated_continuous <- OUwie.sim(
    full_history,
    full_tip_data,
    simmap.tree = TRUE,
    alpha = rep(2.5, 4L),
    sigma.sq = rep(0.7, 4L),
    theta0 = 3.1,
    theta = c(2.5, 2.5, 3.7, 3.7)
  )

  species <- as.character(simulated_continuous[, 1L])
  observed_data <- data.frame(
    species = species,
    state = substr(full_history$states[species], 1L, 1L),
    trait = as.numeric(simulated_continuous[, 2L]),
    stringsAsFactors = FALSE
  )

  list(
    phy = phy,
    data = observed_data,
    full_history = full_history,
    q_matrix = q_matrix,
    truth = c(
      q_observed_A = 0.8,
      q_observed_B = 1.5,
      q_B_to_A = 0.2,
      q_A_to_B = 0.3,
      alpha = 2.5,
      sigma2 = 0.7,
      theta_A = 2.5,
      theta_B = 3.7
    ),
    history_changes = sum(lengths(full_history$maps) - 1L)
  )
}

prepare_problem <- function(simulated) {
  phy <- ape::reorder.phylo(simulated$phy, "pruningwise")
  data <- simulated$data
  rate_cat <- 2L
  root_p <- "yang"
  houwie_data <- organizeHOUwieDat(data, tip.fog = "none", collapse = TRUE)
  n_states <- max(as.numeric(houwie_data$data.cor[, 2L]))
  tree_height <- max(ape::branching.times(phy))
  time_slice <- tree_height + 1
  all_paths <- lapply(
    seq_len(ape::Nnode(phy) + ape::Ntip(phy)),
    function(node) getPathToRoot(phy, node)
  )

  discrete_index <- getDiscreteModel(
    data[, 1:2],
    model = "ER",
    rate.cat = rate_cat,
    dual = FALSE,
    collapse = TRUE
  )
  discrete_index[discrete_index == 0] <- NA
  continuous_index <- getOUParamStructure(
    model = "OUM",
    nObsState = n_states,
    rate.cat = rate_cat,
    null.model = TRUE
  )

  n_transitions <- max(discrete_index, na.rm = TRUE)
  n_alpha <- length(unique(stats::na.omit(continuous_index[1L, ])))
  n_sigma <- length(unique(stats::na.omit(continuous_index[2L, ])))
  n_theta <- length(unique(stats::na.omit(continuous_index[3L, ])))

  lower_original <- c(
    rep(1 / (tree_height * 10000), n_transitions),
    rep(1e-10, n_alpha),
    rep(1e-10, n_sigma),
    rep(min(houwie_data$data.ou[, 3L]) / 10, n_theta)
  )
  upper_alpha <- log(2) / (0.01 * tree_height)
  upper_original <- c(
    rep(1 / (tree_height * 0.0001), n_transitions),
    rep(upper_alpha, n_alpha),
    rep(upper_alpha, n_sigma),
    rep(max(houwie_data$data.ou[, 3L]) * 10, n_theta)
  )

  trait_bin <- cut(houwie_data$data.ou[, 3L], rate_cat, labels = FALSE)
  combinations <- expand.grid(
    seq_len(max(houwie_data$data.cor[, 2L])),
    seq_len(rate_cat)
  )
  discrete_tips <- numeric(length(phy$tip.label))
  for (i in seq_len(nrow(combinations))) {
    matches <- houwie_data$data.cor[, 2L] == combinations[i, 1L] &
      trait_bin == combinations[i, 2L]
    discrete_tips[matches] <- i
  }
  start <- c(
    rep(10 / sum(phy$edge.length), n_transitions),
    rep(log(2) / tree_height, n_alpha),
    rep(log(2) / tree_height, n_sigma),
    getIP.theta(houwie_data$data.ou[, 3L], discrete_tips, continuous_index[3L, ])
  )
  start <- checkStartingUBLB(start, lower_original, upper_original)

  edge_likelihoods <- getEdgeLiks(
    phy,
    houwie_data$data.cor,
    n_states,
    rate_cat,
    time_slice
  )

  raw_objective <- function(log_parameters) {
    value <- tryCatch(
      suppressWarnings(hOUwie.dev(
        p = log_parameters,
        phy = phy,
        data = houwie_data$data.ou,
        rate.cat = rate_cat,
        tip.fog = "none",
        index.disc = discrete_index,
        index.cont = continuous_index,
        root.p = root_p,
        edge_liks_list = edge_likelihoods,
        nSim = n_sim,
        all.paths = all_paths,
        sample_tips = FALSE,
        sample_nodes = FALSE,
        adaptive_sampling = FALSE,
        split.liks = FALSE,
        global_liks_mat = NULL,
        diagn_msg = FALSE
      )),
      error = function(error) 1e10
    )
    if (length(value) != 1L || !is.finite(value)) 1e10 else value
  }

  list(
    raw_objective = raw_objective,
    start = start,
    lower = log(lower_original),
    upper = log(upper_original),
    lower_original = lower_original,
    upper_original = upper_original
  )
}

make_counted_crn_objective <- function(raw_objective, seed) {
  state <- new.env(parent = emptyenv())
  state$evaluations <- 0L
  objective <- makeCommonRandomObjective(
    function(parameters) {
      state$evaluations <- state$evaluations + 1L
      raw_objective(parameters)
    },
    seed
  )
  list(objective = objective, state = state)
}

run_nlopt <- function(algorithm, x0, problem, crn_seed) {
  counted <- make_counted_crn_objective(problem$raw_objective, crn_seed)
  started <- proc.time()[["elapsed"]]
  fit <- nloptr::nloptr(
    x0 = x0,
    eval_f = counted$objective,
    lb = problem$lower,
    ub = problem$upper,
    opts = list(
      algorithm = algorithm,
      maxeval = 0,
      ftol_rel = tolerance,
      xtol_rel = tolerance,
      print_level = 0
    )
  )
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    solution = fit$solution,
    objective = fit$objective,
    evaluations = counted$state$evaluations,
    elapsed = elapsed,
    status = fit$status,
    message = fit$message
  )
}

run_nomad <- function(x0, problem, crn_seed, optimizer_seed) {
  counted <- make_counted_crn_objective(problem$raw_objective, crn_seed)
  started <- proc.time()[["elapsed"]]
  fit <- crs::snomadr(
    eval.f = counted$objective,
    n = length(x0),
    bbin = rep(0L, length(x0)),
    bbout = 0L,
    x0 = x0,
    lb = problem$lower,
    ub = problem$upper,
    nmulti = 0L,
    random.seed = optimizer_seed,
    opts = list(
      "MAX_BB_EVAL" = 1000000L,
      "INITIAL_MESH_SIZE" = "r1.0e-1",
      "MIN_MESH_SIZE" = paste0("r", format(tolerance, scientific = TRUE)),
      "DISPLAY_DEGREE" = 0L
    ),
    display.nomad.progress = FALSE,
    snomadr.environment = environment(counted$objective)
  )
  elapsed <- proc.time()[["elapsed"]] - started

  solution <- fit$solution
  if (is.null(solution)) solution <- fit$xbest
  if (is.null(solution)) solution <- fit$x
  objective <- fit$objective
  if (is.null(objective)) objective <- fit$fbest
  if (is.null(objective)) objective <- counted$objective(solution)

  list(
    solution = as.numeric(solution),
    objective = as.numeric(objective)[1L],
    evaluations = counted$state$evaluations,
    elapsed = elapsed,
    status = if (is.null(fit$status)) NA_integer_ else fit$status,
    message = if (is.null(fit$message)) "" else as.character(fit$message),
    nomad_result_names = names(fit)
  )
}

simulated <- simulate_challenging_data(n_tips = n_tips)
problem <- prepare_problem(simulated)
set.seed(benchmark_seed + 2L)
starts <- generateMultiStarting(
  problem$start,
  getDiscreteModel(simulated$data[, 1:2], "ER", 2L, FALSE, TRUE),
  getOUParamStructure("OUM", 2L, 2L, TRUE),
  n_starts,
  problem$lower_original,
  problem$upper_original
)
starts <- lapply(starts, log)
crn_seeds <- sample.int(.Machine$integer.max, n_starts)
optimizer_seeds <- sample.int(.Machine$integer.max, n_starts)
validation_seeds <- sample.int(.Machine$integer.max, n_validation)

saveRDS(
  list(
    simulated = simulated,
    starts = starts,
    crn_seeds = crn_seeds,
    optimizer_seeds = optimizer_seeds,
    validation_seeds = validation_seeds,
    n_sim = n_sim,
    n_tips = n_tips,
    optimizers = optimizer_names,
    tolerance = tolerance
  ),
  file.path(results_dir, "benchmark-input.rds")
)

tasks <- expand.grid(
  optimizer = optimizer_names,
  start_id = seq_len(n_starts),
  stringsAsFactors = FALSE
)

run_task <- function(task_id) {
  task <- tasks[task_id, ]
  output_file <- file.path(
    fits_dir,
    paste0(tolower(task$optimizer), "-start-", task$start_id, ".rds")
  )
  if (file.exists(output_file)) return(readRDS(output_file))

  x0 <- starts[[task$start_id]]
  if (task$optimizer == "SBPLX") {
    fit <- run_nlopt("NLOPT_LN_SBPLX", x0, problem, crn_seeds[task$start_id])
  } else if (task$optimizer == "COBYLA") {
    fit <- run_nlopt("NLOPT_LN_COBYLA", x0, problem, crn_seeds[task$start_id])
  } else {
    fit <- run_nomad(
      x0,
      problem,
      crn_seeds[task$start_id],
      optimizer_seeds[task$start_id]
    )
  }

  fit$optimizer <- task$optimizer
  fit$start_id <- task$start_id
  fit$crn_seed <- crn_seeds[task$start_id]
  fit$training_loglik <- -fit$objective
  saveRDS(fit, output_file)
  fit
}

cores <- suppressWarnings(as.integer(Sys.getenv("HOUWIE_BENCH_CORES", "3")))
if (!is.finite(cores) || cores < 1L) cores <- 1L
fits <- parallel::mclapply(
  seq_len(nrow(tasks)),
  run_task,
  mc.cores = min(cores, nrow(tasks)),
  mc.preschedule = FALSE
)

if (any(vapply(fits, inherits, logical(1), what = "try-error"))) {
  stop("At least one optimizer task failed; completed task files can be resumed.")
}

validate_fit <- function(fit) {
  validation_loglik <- vapply(
    validation_seeds,
    function(seed) {
      -withHOUwieSeed(seed, problem$raw_objective(fit$solution))
    },
    numeric(1)
  )
  fit$validation_loglik <- validation_loglik
  fit$validation_mean <- mean(validation_loglik)
  fit$validation_median <- stats::median(validation_loglik)
  fit$validation_sd <- stats::sd(validation_loglik)
  fit$optimization_optimism <- fit$training_loglik - fit$validation_mean
  fit$parameters <- exp(fit$solution)
  fit
}

fits <- lapply(fits, validate_fit)
saveRDS(fits, file.path(results_dir, "validated-fits.rds"))

summary_rows <- lapply(fits, function(fit) {
  data.frame(
    optimizer = fit$optimizer,
    start_id = fit$start_id,
    training_loglik = fit$training_loglik,
    validation_mean = fit$validation_mean,
    validation_median = fit$validation_median,
    validation_sd = fit$validation_sd,
    optimization_optimism = fit$optimization_optimism,
    evaluations = fit$evaluations,
    elapsed_seconds = fit$elapsed,
    status = paste(fit$status, collapse = ";"),
    message = paste(fit$message, collapse = ";"),
    stringsAsFactors = FALSE
  )
})
fit_summary <- do.call(rbind, summary_rows)
utils::write.csv(
  fit_summary,
  file.path(results_dir, "fit-summary.csv"),
  row.names = FALSE
)

parameter_names <- names(simulated$truth)
parameter_rows <- lapply(fits, function(fit) {
  values <- fit$parameters
  names(values) <- parameter_names
  data.frame(
    optimizer = fit$optimizer,
    start_id = fit$start_id,
    parameter = parameter_names,
    estimate = values,
    truth = simulated$truth,
    stringsAsFactors = FALSE
  )
})
utils::write.csv(
  do.call(rbind, parameter_rows),
  file.path(results_dir, "parameter-estimates.csv"),
  row.names = FALSE
)

aggregate_summary <- do.call(
  rbind,
  lapply(split(fit_summary, fit_summary$optimizer), function(rows) {
    data.frame(
      optimizer = rows$optimizer[1L],
      mean_validation_loglik = mean(rows$validation_mean),
      median_validation_loglik = stats::median(rows$validation_mean),
      between_start_sd = stats::sd(rows$validation_mean),
      mean_evaluations = mean(rows$evaluations),
      mean_elapsed_seconds = mean(rows$elapsed_seconds),
      mean_optimism = mean(rows$optimization_optimism),
      stringsAsFactors = FALSE
    )
  })
)
utils::write.csv(
  aggregate_summary,
  file.path(results_dir, "aggregate-summary.csv"),
  row.names = FALSE
)

print(aggregate_summary, row.names = FALSE)
