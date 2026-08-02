#!/usr/bin/env Rscript

# Prototype a proposal-corrected likelihood for hOUwie. Histories are sampled
# with replacement from either the discrete-only or continuous-informed node
# proposal. Each contribution is divided by its exact sequential proposal
# probability before averaging.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]),
                             mustWork = TRUE)
benchmark_dir <- dirname(script_file)
repo_root <- normalizePath(file.path(benchmark_dir, "..", "..", ".."),
                           mustWork = TRUE)

if(!requireNamespace("pkgload", quietly = TRUE)){
  stop("pkgload is required to load the current local OUwie source.",
       call. = FALSE)
}
pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                  helpers = FALSE)

fixture <- readRDS(file.path(
  benchmark_dir, "results", "scaling", "fixture-tips300-seed8.rds"
))

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

log_sum_exp <- function(x){
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

proposal_log_probability <- function(state_sample, Pj, root_state,
                                     tree_plan){
  root_sample <- state_sample[[tree_plan$root.edges[1]]][1]
  value <- log(root_state[root_sample] / sum(root_state))
  for(edge_i in tree_plan$rev.pruning.order){
    from <- state_sample[[edge_i]][1]
    count <- 2L
    n_nodes <- length(state_sample[[edge_i]])
    for(interval_i in (n_nodes - 1L):1L){
      to <- state_sample[[edge_i]][count]
      probabilities <- Pj[[edge_i]][from, , interval_i]
      value <- value + log(probabilities[to] / sum(probabilities))
      from <- to
      count <- count + 1L
    }
  }
  value
}

phy <- ape::reorder.phylo(fixture$phy, "pruningwise")
data <- OUwie:::matchTipsAndData(phy$tip.label, fixture$data)
houwie_data <- OUwie:::organizeHOUwieDat(data, "none", TRUE)
edge_liks_initial <- OUwie:::getEdgeLiks(
  phy,
  houwie_data$data.cor,
  2,
  1,
  max(ape::branching.times(phy)) + 1
)
tree_plan <- OUwie:::getHOUwieTreePlan(phy, edge_liks_initial)

prepare_proposal <- function(candidate, continuous_informed){
  Q <- matrix(
    c(-candidate$q, candidate$q, candidate$q, -candidate$q),
    2, 2, byrow = TRUE
  )
  rate_matrix <- rbind(
    alpha = rep(candidate$alpha, 2),
    sigma.sq = rep(candidate$sigma_sq, 2),
    theta = unname(c(candidate$theta_a, candidate$theta_b))
  )
  edge_liks <- edge_liks_initial
  if(continuous_informed){
    edge_liks <- OUwie:::getCherryConditionals(
      phy,
      houwie_data$data.ou,
      rate_matrix,
      Q,
      edge_liks_initial,
      tree_plan$tip.paths
    )
  }
  for(edge_i in seq_along(edge_liks)){
    edge_liks[[edge_i]] <- edge_liks[[edge_i]] *
      edge_liks_initial[[edge_i]]
  }
  conditional <- OUwie:::getConditionalInternodeLik(phy, Q, edge_liks)
  root_prior <- c(0.5, 0.5)
  sampler <- OUwie:::getInternodeMap(
    phy,
    Q,
    conditional$edge_liks_list,
    conditional$root_state,
    root_prior,
    nSim = 1,
    check_vector = NA,
    max.attempts = 2,
    tree.plan = tree_plan
  )
  list(
    Q = Q,
    rate_matrix = rate_matrix,
    root_prior = root_prior,
    root_state = conditional$root_state,
    Pij = sampler$Pij,
    Pj = sampler$Pj
  )
}

draw_history <- function(proposal){
  OUwie:::getInternodeStateSample(
    proposal$Pj,
    proposal$root_state,
    tree_plan$root.edges,
    tree_plan$rev.pruning.order,
    tree_plan$edge.index,
    2,
    tree_plan$number.of.nodes.per.edge,
    tree_plan$child.edges
  )
}

evaluate_importance_sample <- function(candidate, continuous_informed, n_sim){
  proposal <- prepare_proposal(candidate, continuous_informed)
  continuous_loglik <- target_discrete_loglik <- proposal_loglik <-
    numeric(n_sim)
  history_ids <- character(n_sim)
  for(draw_i in seq_len(n_sim)){
    state_sample <- draw_history(proposal)
    compact_map <- OUwie:::getMapFromStateSample(
      tree_plan$map.template, state_sample
    )
    continuous_loglik[draw_i] <- OUwie:::OUwie.basic(
      phy,
      houwie_data$data.ou,
      simmap.tree = TRUE,
      scaleHeight = FALSE,
      alpha = proposal$rate_matrix["alpha", ],
      sigma.sq = proposal$rate_matrix["sigma.sq", ],
      theta = proposal$rate_matrix["theta", ],
      algorithm = "three.point",
      tip.paths = tree_plan$tip.paths,
      tip.fog = "none",
      map = compact_map,
      tree.plan = tree_plan
    )
    target_discrete_loglik[draw_i] <- OUwie:::getStateSampleProb(
      state_sample,
      proposal$Pij,
      proposal$root_prior,
      tree_plan$root.edges
    )
    proposal_loglik[draw_i] <- proposal_log_probability(
      state_sample,
      proposal$Pj,
      proposal$root_state,
      tree_plan
    )
    history_ids[draw_i] <- OUwie:::stateSampleId(state_sample)
  }
  log_weights <- target_discrete_loglik + continuous_loglik -
    proposal_loglik
  normalized_weights <- exp(log_weights - log_sum_exp(log_weights))
  unique_draw <- !duplicated(history_ids)
  list(
    importance_loglik = log_sum_exp(log_weights) - log(n_sim),
    uncorrected_unique_logsum = log_sum_exp(
      target_discrete_loglik[unique_draw] + continuous_loglik[unique_draw]
    ),
    effective_sample_size = 1 / sum(normalized_weights^2),
    maximum_weight_fraction = max(normalized_weights),
    unique_histories = sum(unique_draw),
    discrete_minus_proposal_sd = sd(
      target_discrete_loglik - proposal_loglik
    )
  )
}

repetitions <- suppressWarnings(as.integer(
  Sys.getenv("HOUWIE_IMPORTANCE_REPS", "3")
))
counts <- as.integer(strsplit(
  Sys.getenv("HOUWIE_IMPORTANCE_NSIM", "100,500"), ",", fixed = TRUE
)[[1]])
if(is.na(repetitions) || repetitions < 1L ||
   any(is.na(counts)) || any(counts < 1L)){
  stop("Importance-sampling repetitions and counts must be positive integers.",
       call. = FALSE)
}

results_dir <- file.path(
  benchmark_dir, "results", "history-integration-diagnosis"
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(results_dir, "importance-sampling.csv")
if(identical(Sys.getenv("HOUWIE_BENCH_OVERWRITE"), "1")) unlink(results_path)

read_results <- function(){
  if(file.exists(results_path)) {
    read.csv(results_path, stringsAsFactors = FALSE)
  } else NULL
}

for(continuous_informed in c(FALSE, TRUE)){
  proposal_name <- if(continuous_informed) {
    "continuous-informed-nodes"
  } else "discrete-only"
  for(n_sim in counts){
    for(repetition in seq_len(repetitions)){
      run_seed <- 99000L + 10000L * as.integer(continuous_informed) +
        1000L * n_sim + repetition
      for(candidate_i in seq_len(nrow(candidates))){
        candidate <- candidates[candidate_i, ]
        candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidate$candidate))
        task_id <- paste("importance", proposal_name, n_sim, candidate_id,
                         repetition, sep = "-")
        old <- read_results()
        if(!is.null(old) && task_id %in% old$task_id) next
        cat("Running", task_id, "\n")
        set.seed(run_seed)
        error_message <- NA_character_
        started <- proc.time()[["elapsed"]]
        result <- tryCatch(
          evaluate_importance_sample(
            candidate, continuous_informed, n_sim
          ),
          error = function(e){
            error_message <<- conditionMessage(e)
            NULL
          }
        )
        elapsed <- proc.time()[["elapsed"]] - started
        row <- data.frame(
          task_id = task_id,
          proposal = proposal_name,
          n_sim = n_sim,
          repetition = repetition,
          seed = run_seed,
          candidate = candidate$candidate,
          importance_loglik = if(is.null(result)) NA_real_ else
            result$importance_loglik,
          uncorrected_unique_logsum = if(is.null(result)) NA_real_ else
            result$uncorrected_unique_logsum,
          effective_sample_size = if(is.null(result)) NA_real_ else
            result$effective_sample_size,
          maximum_weight_fraction = if(is.null(result)) NA_real_ else
            result$maximum_weight_fraction,
          unique_histories = if(is.null(result)) NA_integer_ else
            result$unique_histories,
          discrete_minus_proposal_sd = if(is.null(result)) NA_real_ else
            result$discrete_minus_proposal_sd,
          elapsed_seconds = elapsed,
          error = error_message,
          stringsAsFactors = FALSE
        )
        combined <- if(is.null(old)) row else rbind(old, row)
        write.csv(combined, results_path, row.names = FALSE)
      }
    }
  }
}

cat("Finished importance-sampling prototype\n")
cat("Results:", results_path, "\n")
