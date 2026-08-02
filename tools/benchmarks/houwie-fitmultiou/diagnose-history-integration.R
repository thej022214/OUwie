#!/usr/bin/env Rscript

# Diagnose why hOUwie and fitmultiOU rank the same parameter vectors
# differently on Liam Revell's 300-tip simulation. The probes vary one part of
# the history integration at a time and save after every evaluation so a run
# can be resumed.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]),
                             mustWork = TRUE)
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

if(!requireNamespace("pkgload", quietly = TRUE)){
  stop("pkgload is required to load the current local OUwie source.",
       call. = FALSE)
}
if(!requireNamespace("phytools", quietly = TRUE)){
  stop("phytools is required for the discrete-trait reference likelihood.",
       call. = FALSE)
}
pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                  helpers = FALSE)

fixture_path <- file.path(
  benchmark_dir, "results", "scaling", "fixture-tips300-seed8.rds"
)
if(!file.exists(fixture_path)){
  stop("The exact 300-tip fixture is missing. Run the scaling profile first.",
       call. = FALSE)
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

candidate_parameters <- function(candidate){
  unname(c(
    candidate$q,
    candidate$alpha,
    candidate$sigma_sq,
    candidate$theta_a,
    candidate$theta_b
  ))
}

results_dir <- file.path(
  benchmark_dir, "results", "history-integration-diagnosis"
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

repetitions <- suppressWarnings(as.integer(
  Sys.getenv("HOUWIE_DIAG_REPS", "3")
))
if(is.na(repetitions) || repetitions < 1L){
  stop("HOUWIE_DIAG_REPS must be a positive integer.", call. = FALSE)
}

overwrite <- identical(Sys.getenv("HOUWIE_BENCH_OVERWRITE"), "1")

append_task <- function(path, row){
  old <- if(file.exists(path)) {
    read.csv(path, stringsAsFactors = FALSE)
  } else NULL
  if(!is.null(old) && row$task_id %in% old$task_id) return(invisible(FALSE))
  combined <- if(is.null(old)) row else rbind(old, row)
  write.csv(combined, path, row.names = FALSE)
  invisible(TRUE)
}

task_done <- function(path, task_id){
  if(!file.exists(path)) return(FALSE)
  task_id %in% read.csv(path, stringsAsFactors = FALSE)$task_id
}

discrete_loglik <- local({
  cache <- new.env(parent = emptyenv())
  function(q){
    key <- format(q, digits = 16)
    if(exists(key, cache, inherits = FALSE)) return(get(key, cache))
    Q <- matrix(
      c(-q, q, q, -q), 2, 2, byrow = TRUE,
      dimnames = list(c("a", "b"), c("a", "b"))
    )
    value <- as.numeric(phytools::fitMk(
      fixture$phy,
      fixture$discrete,
      fixedQ = Q,
      pi = c(a = 0.5, b = 0.5)
    )$logLik)
    assign(key, value, cache)
    value
  }
})

result_row <- function(task_id, probe, setting, candidate, repetition,
                       seed, resolution, elapsed, fit, error_message){
  n_maps <- if(is.null(fit)) NA_integer_ else length(fit$all_cont_liks)
  continuous_loglik <- if(is.null(fit)) NA_real_ else fit$ContLik
  corrected_mc_loglik <- if(is.null(fit) || !n_maps) NA_real_ else
    discrete_loglik(candidate$q) + continuous_loglik - log(n_maps)
  data.frame(
    task_id = task_id,
    probe = probe,
    setting = setting,
    candidate = candidate$candidate,
    repetition = repetition,
    seed = seed,
    resolution = resolution,
    retained_histories = n_maps,
    elapsed_seconds = elapsed,
    current_loglik = if(is.null(fit)) NA_real_ else fit$loglik,
    discrete_tip_loglik = if(is.null(fit)) NA_real_ else
      discrete_loglik(candidate$q),
    corrected_mc_loglik = corrected_mc_loglik,
    sampled_discrete_logsum = if(is.null(fit)) NA_real_ else fit$DiscLik,
    sampled_continuous_logsum = continuous_loglik,
    error = error_message,
    stringsAsFactors = FALSE
  )
}

run_public_probe <- function(path, probe, settings){
  if(overwrite && file.exists(path)) unlink(path)
  for(setting_i in seq_len(nrow(settings))){
    setting <- settings[setting_i, ]
    for(repetition in seq_len(repetitions)){
      run_seed <- 94000L + 100L * setting_i + repetition
      for(candidate_i in seq_len(nrow(candidates))){
        candidate <- candidates[candidate_i, ]
        candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidate$candidate))
        task_id <- paste(probe, setting$setting, candidate_id, repetition,
                         sep = "-")
        if(task_done(path, task_id)) next
        cat("Running", task_id, "\n")
        set.seed(run_seed)
        error_message <- NA_character_
        started <- proc.time()[["elapsed"]]
        fit <- tryCatch(
          hOUwie(
            phy = fixture$phy,
            data = fixture$data,
            rate.cat = 1,
            discrete_model = "ER",
            continuous_model = "OUM",
            nSim = setting$n_sim,
            p = candidate_parameters(candidate),
            sample_tips = FALSE,
            sample_nodes = setting$sample_nodes,
            adaptive_sampling = setting$adaptive_sampling,
            common_random_numbers = TRUE,
            quiet = TRUE
          ),
          error = function(e){
            error_message <<- conditionMessage(e)
            NULL
          }
        )
        elapsed <- proc.time()[["elapsed"]] - started
        append_task(path, result_row(
          task_id = task_id,
          probe = probe,
          setting = setting$setting,
          candidate = candidate,
          repetition = repetition,
          seed = run_seed,
          resolution = setting$n_sim,
          elapsed = elapsed,
          fit = fit,
          error_message = error_message
        ))
      }
    }
  }
}

mode_path <- file.path(results_dir, "sampling-modes.csv")
sampling_modes <- data.frame(
  setting = c(
    "discrete-only",
    "continuous-informed-nodes",
    "adaptive",
    "continuous-informed-and-adaptive"
  ),
  n_sim = 100L,
  sample_nodes = c(FALSE, TRUE, FALSE, TRUE),
  adaptive_sampling = c(FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
run_public_probe(mode_path, "sampling-mode", sampling_modes)

count_path <- file.path(results_dir, "history-counts.csv")
history_counts <- data.frame(
  setting = paste0("nSim-", c(25L, 100L, 500L)),
  n_sim = c(25L, 100L, 500L),
  sample_nodes = FALSE,
  adaptive_sampling = FALSE,
  stringsAsFactors = FALSE
)
run_public_probe(count_path, "history-count", history_counts)

# The public function currently places only the two endpoint states on an edge,
# with a change halfway along the edge when they differ. This probe inserts
# progressively more latent states along longer edges while leaving every other
# likelihood component unchanged.
resolution_path <- file.path(results_dir, "within-edge-resolution.csv")
if(overwrite && file.exists(resolution_path)) unlink(resolution_path)

phy <- ape::reorder.phylo(fixture$phy, "pruningwise")
data <- OUwie:::matchTipsAndData(phy$tip.label, fixture$data)
houwie_data <- OUwie:::organizeHOUwieDat(data, "none", TRUE)
index_disc <- OUwie:::getDiscreteModel(
  houwie_data$data.cor, "ER", 1, FALSE, TRUE
)
index_disc[index_disc == 0] <- NA
index_cont <- OUwie:::getOUParamStructure("OUM", 2, 1, FALSE)
tree_height <- max(ape::branching.times(phy))
slice_settings <- data.frame(
  setting = c("one-midpoint-per-edge", "two-time-unit-grid",
              "one-time-unit-grid", "half-time-unit-grid"),
  time_slice = c(tree_height + 1, 2, 1, 0.5),
  stringsAsFactors = FALSE
)

for(setting_i in seq_len(nrow(slice_settings))){
  setting <- slice_settings[setting_i, ]
  edge_liks <- OUwie:::getEdgeLiks(
    phy, houwie_data$data.cor, 2, 1, setting$time_slice
  )
  tree_plan <- OUwie:::getHOUwieTreePlan(phy, edge_liks)
  for(repetition in seq_len(repetitions)){
    run_seed <- 95000L + 100L * setting_i + repetition
    for(candidate_i in seq_len(nrow(candidates))){
      candidate <- candidates[candidate_i, ]
      candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidate$candidate))
      task_id <- paste("within-edge", setting$setting, candidate_id,
                       repetition, sep = "-")
      if(task_done(resolution_path, task_id)) next
      cat("Running", task_id, "\n")
      set.seed(run_seed)
      error_message <- NA_character_
      started <- proc.time()[["elapsed"]]
      raw_fit <- tryCatch(
        OUwie:::hOUwie.dev(
          p = log(candidate_parameters(candidate)),
          phy = phy,
          data = houwie_data$data.ou,
          rate.cat = 1,
          tip.fog = "none",
          index.disc = index_disc,
          index.cont = index_cont,
          root.p = "yang",
          edge_liks_list = edge_liks,
          nSim = 100,
          all.paths = tree_plan$all.paths,
          sample_tips = FALSE,
          sample_nodes = FALSE,
          adaptive_sampling = FALSE,
          split.liks = TRUE,
          global_liks_mat = NULL,
          diagn_msg = FALSE,
          tree.plan = tree_plan
        ),
        error = function(e){
          error_message <<- conditionMessage(e)
          NULL
        }
      )
      elapsed <- proc.time()[["elapsed"]] - started
      fit <- if(is.null(raw_fit)) NULL else list(
        loglik = raw_fit$TotalLik,
        DiscLik = raw_fit$DiscLik,
        ContLik = raw_fit$ContLik,
        all_cont_liks = raw_fit$llik_continuous
      )
      row <- result_row(
        task_id = task_id,
        probe = "within-edge-resolution",
        setting = setting$setting,
        candidate = candidate,
        repetition = repetition,
        seed = run_seed,
        resolution = setting$time_slice,
        elapsed = elapsed,
        fit = fit,
        error_message = error_message
      )
      row$mean_latent_intervals_per_edge <- mean(
        tree_plan$number.of.edges.per.edge
      )
      append_task(resolution_path, row)
    }
  }
}

# On a six-tip tree there are only 2^5 = 32 assignments of states to internal
# nodes. Enumerating all of them gives the exact likelihood for hOUwie's
# midpoint-history model and lets us measure the error from retaining only a
# sampled subset of histories.
small_exact_path <- file.path(results_dir, "small-tree-exact-enumeration.csv")
small_sampling_path <- file.path(results_dir, "small-tree-sampled-subsets.csv")
if(overwrite){
  unlink(c(small_exact_path, small_sampling_path))
}

set.seed(96001)
small_phy <- phytools::pbtree(n = 6, scale = 10)
truth_Q <- matrix(c(-0.3, 0.3, 0.3, -0.3), 2, 2, byrow = TRUE)
small_simulation <- hOUwie.sim(
  small_phy,
  Q = truth_Q,
  root.freqs = c(1, 0),
  alpha = c(0.5, 0.5),
  sigma.sq = c(0.2, 0.2),
  theta0 = 5,
  theta = c(5, 10)
)
small_data <- small_simulation$data
small_phy <- ape::reorder.phylo(small_phy, "pruningwise")
small_data <- OUwie:::matchTipsAndData(small_phy$tip.label, small_data)
small_houwie_data <- OUwie:::organizeHOUwieDat(small_data, "none", TRUE)
small_edge_liks <- OUwie:::getEdgeLiks(
  small_phy,
  small_houwie_data$data.cor,
  2,
  1,
  max(ape::branching.times(small_phy)) + 1
)
small_tree_plan <- OUwie:::getHOUwieTreePlan(small_phy, small_edge_liks)
small_tip_states <- as.integer(small_houwie_data$data.cor[, 2])
names(small_tip_states) <- small_houwie_data$data.cor[, 1]
small_tip_states <- small_tip_states[small_phy$tip.label]
small_internal_nodes <- (ape::Ntip(small_phy) + 1):(
  ape::Ntip(small_phy) + ape::Nnode(small_phy)
)
state_assignments <- expand.grid(rep(list(1:2), length(small_internal_nodes)))

log_sum_exp <- function(x){
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

exact_midpoint_loglik <- function(candidate){
  q <- candidate$q
  Q <- matrix(c(-q, q, q, -q), 2, 2, byrow = TRUE)
  edge_transitions <- lapply(
    small_phy$edge.length,
    function(edge_length) expm::expm(Q * edge_length)
  )
  joint_loglik <- numeric(nrow(state_assignments))
  for(assignment_i in seq_len(nrow(state_assignments))){
    states <- integer(ape::Ntip(small_phy) + ape::Nnode(small_phy))
    states[seq_len(ape::Ntip(small_phy))] <- small_tip_states
    states[small_internal_nodes] <- as.integer(
      state_assignments[assignment_i, ]
    )
    discrete_part <- log(0.5) + sum(vapply(
      seq_len(nrow(small_phy$edge)),
      function(edge_i){
        log(edge_transitions[[edge_i]][
          states[small_phy$edge[edge_i, 1]],
          states[small_phy$edge[edge_i, 2]]
        ])
      },
      numeric(1)
    ))
    midpoint_map <- mapply(
      function(ancestor, descendant, edge_length){
        setNames(
          rep(edge_length / 2, 2),
          c(states[ancestor], states[descendant])
        )
      },
      ancestor = small_phy$edge[, 1],
      descendant = small_phy$edge[, 2],
      edge_length = small_phy$edge.length,
      SIMPLIFY = FALSE
    )
    continuous_part <- OUwie:::OUwie.basic(
      small_phy,
      small_houwie_data$data.ou,
      simmap.tree = TRUE,
      scaleHeight = FALSE,
      alpha = rep(candidate$alpha, 2),
      sigma.sq = rep(candidate$sigma_sq, 2),
      theta = unname(c(candidate$theta_a, candidate$theta_b)),
      algorithm = "three.point",
      tip.paths = small_tree_plan$tip.paths,
      tip.fog = "none",
      map = midpoint_map,
      tree.plan = small_tree_plan
    )
    joint_loglik[assignment_i] <- discrete_part + continuous_part
  }
  log_sum_exp(joint_loglik)
}

if(!file.exists(small_exact_path)){
  exact_rows <- lapply(seq_len(nrow(candidates)), function(candidate_i){
    candidate <- candidates[candidate_i, ]
    data.frame(
      candidate = candidate$candidate,
      possible_histories = nrow(state_assignments),
      exact_midpoint_loglik = exact_midpoint_loglik(candidate),
      stringsAsFactors = FALSE
    )
  })
  write.csv(do.call(rbind, exact_rows), small_exact_path, row.names = FALSE)
}
small_exact <- read.csv(small_exact_path, stringsAsFactors = FALSE)

small_counts <- c(5L, 15L, 100L)
small_repetitions <- max(30L, repetitions)
for(n_sim in small_counts){
  for(repetition in seq_len(small_repetitions)){
    run_seed <- 97000L + 1000L * n_sim + repetition
    for(candidate_i in seq_len(nrow(candidates))){
      candidate <- candidates[candidate_i, ]
      candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidate$candidate))
      task_id <- paste("small-tree", n_sim, candidate_id, repetition,
                       sep = "-")
      if(task_done(small_sampling_path, task_id)) next
      set.seed(run_seed)
      error_message <- NA_character_
      fit <- tryCatch(
        hOUwie(
          phy = small_phy,
          data = small_data,
          rate.cat = 1,
          discrete_model = "ER",
          continuous_model = "OUM",
          nSim = n_sim,
          p = candidate_parameters(candidate),
          sample_tips = FALSE,
          sample_nodes = FALSE,
          adaptive_sampling = FALSE,
          common_random_numbers = FALSE,
          quiet = TRUE
        ),
        error = function(e){
          error_message <<- conditionMessage(e)
          NULL
        }
      )
      exact_value <- small_exact$exact_midpoint_loglik[
        small_exact$candidate == candidate$candidate
      ]
      row <- data.frame(
        task_id = task_id,
        n_sim = n_sim,
        repetition = repetition,
        seed = run_seed,
        candidate = candidate$candidate,
        retained_histories = if(is.null(fit)) NA_integer_ else
          length(fit$all_cont_liks),
        sampled_loglik = if(is.null(fit)) NA_real_ else fit$loglik,
        exact_midpoint_loglik = exact_value,
        error_from_exact = if(is.null(fit)) NA_real_ else
          fit$loglik - exact_value,
        error = error_message,
        stringsAsFactors = FALSE
      )
      append_task(small_sampling_path, row)
    }
  }
}

# Changing the history proposal must not change the target likelihood. On the
# small tree the exact target above is fixed, so any proposal-dependent shift is
# approximation error rather than biological signal.
small_proposal_path <- file.path(
  results_dir, "small-tree-proposal-dependence.csv"
)
if(overwrite && file.exists(small_proposal_path)) unlink(small_proposal_path)
small_proposals <- data.frame(
  proposal = c("discrete-only", "continuous-informed-nodes"),
  sample_nodes = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
for(proposal_i in seq_len(nrow(small_proposals))){
  proposal <- small_proposals[proposal_i, ]
  for(n_sim in c(5L, 15L)){
    for(repetition in seq_len(small_repetitions)){
      run_seed <- 98000L + 10000L * proposal_i + 1000L * n_sim + repetition
      for(candidate_i in seq_len(nrow(candidates))){
        candidate <- candidates[candidate_i, ]
        candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidate$candidate))
        task_id <- paste("small-proposal", proposal$proposal, n_sim,
                         candidate_id, repetition, sep = "-")
        if(task_done(small_proposal_path, task_id)) next
        set.seed(run_seed)
        error_message <- NA_character_
        fit <- tryCatch(
          hOUwie(
            phy = small_phy,
            data = small_data,
            rate.cat = 1,
            discrete_model = "ER",
            continuous_model = "OUM",
            nSim = n_sim,
            p = candidate_parameters(candidate),
            sample_tips = FALSE,
            sample_nodes = proposal$sample_nodes,
            adaptive_sampling = FALSE,
            common_random_numbers = FALSE,
            quiet = TRUE
          ),
          error = function(e){
            error_message <<- conditionMessage(e)
            NULL
          }
        )
        exact_value <- small_exact$exact_midpoint_loglik[
          small_exact$candidate == candidate$candidate
        ]
        row <- data.frame(
          task_id = task_id,
          proposal = proposal$proposal,
          n_sim = n_sim,
          repetition = repetition,
          seed = run_seed,
          candidate = candidate$candidate,
          retained_histories = if(is.null(fit)) NA_integer_ else
            length(fit$all_cont_liks),
          sampled_loglik = if(is.null(fit)) NA_real_ else fit$loglik,
          exact_midpoint_loglik = exact_value,
          error_from_exact = if(is.null(fit)) NA_real_ else
            fit$loglik - exact_value,
          error = error_message,
          stringsAsFactors = FALSE
        )
        append_task(small_proposal_path, row)
      }
    }
  }
}

# fitmultiOU offers several root treatments. Re-evaluating the candidates under
# each one checks whether the cross-method ranking is a root-prior artifact.
root_path <- file.path(results_dir, "fitmultiou-root-treatments.csv")
if(overwrite && file.exists(root_path)) unlink(root_path)
if("fitmultiOU" %in% getNamespaceExports("phytools")){
  initial <- setNames(
    c(5, 10, 0.5, 0.2, 0.3),
    c("theta[a]", "theta[b]", "alpha", "sigsq", "q[1]")
  )
  for(root in c("mle", "flat", "nuisance")){
    setup_id <- paste("fitmulti-root", root, sep = "-")
    missing_candidates <- vapply(seq_len(nrow(candidates)), function(i){
      candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidates$candidate[i]))
      !task_done(root_path, paste(setup_id, candidate_id, sep = "-"))
    }, logical(1))
    if(!any(missing_candidates)) next
    fitmulti_setup <- phytools::fitmultiOU(
      fixture$phy,
      fixture$trait,
      fixture$discrete,
      model = "ER",
      levs = 100,
      parallel = FALSE,
      root = root,
      trace = 0,
      maxit = 0,
      rand_start = FALSE,
      init = initial
    )
    for(candidate_i in which(missing_candidates)){
      candidate <- candidates[candidate_i, ]
      candidate_id <- gsub("[^a-z0-9]+", "-", tolower(candidate$candidate))
      task_id <- paste(setup_id, candidate_id, sep = "-")
      value <- fitmulti_setup$lik(
        theta = unname(c(candidate$theta_a, candidate$theta_b)),
        alpha = log(candidate$alpha),
        sigsq = log(candidate$sigma_sq),
        q = log(candidate$q)
      )
      append_task(root_path, data.frame(
        task_id = task_id,
        root = root,
        candidate = candidate$candidate,
        loglik = value,
        stringsAsFactors = FALSE
      ))
    }
  }
}else{
  cat("Skipping root-treatment probe: fitmultiOU is not installed.\n")
}

cat("Finished history-integration diagnosis\n")
cat("Results:", results_dir, "\n")
