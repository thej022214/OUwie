#!/usr/bin/env Rscript

# Independent continuous-likelihood audit for master and houwie-update. The
# reference calculation propagates linear-Gaussian OU moments directly along
# each painted edge and does not call OUwie's likelihood implementation.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."),
                           mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)

reorder_simmap <- function(simmap, order){
  index <- ape::reorder.phylo(simmap, order, index.only = TRUE)
  simmap$edge <- simmap$edge[index, , drop = FALSE]
  simmap$edge.length <- simmap$edge.length[index]
  simmap$maps <- simmap$maps[index]
  simmap$mapped.edge <- simmap$mapped.edge[index, , drop = FALSE]
  attr(simmap, "order") <- order
  simmap
}

direct_ou_loglik <- function(simmap, traits, alpha, sigma_sq, theta){
  state_names <- names(alpha)
  n_tip <- ape::Ntip(simmap)
  n_all <- n_tip + simmap$Nnode
  root <- n_tip + 1L
  root_edges <- which(simmap$edge[, 1] == root)
  root_states <- vapply(simmap$maps[root_edges], function(x){
    names(x)[1]
  }, character(1))
  if(length(unique(root_states)) != 1L){
    stop("Root-descending maps disagree about the root state")
  }
  root_state <- match(root_states[1], state_names)

  node_mean <- rep(NA_real_, n_all)
  node_cov <- matrix(0, n_all, n_all)
  assigned <- rep(FALSE, n_all)
  node_mean[root] <- theta[root_state]
  # hOUwie's default parameterization fixes the root at its regime optimum.
  node_cov[root, root] <- 0
  assigned[root] <- TRUE

  edge_order <- ape::reorder.phylo(simmap, "cladewise", index.only = TRUE)
  for(edge_i in edge_order){
    ancestor <- simmap$edge[edge_i, 1]
    descendant <- simmap$edge[edge_i, 2]
    if(!assigned[ancestor]){
      stop("Traversal reached a child before its ancestor")
    }
    edge_a <- 1
    edge_b <- 0
    edge_v <- 0
    edge_map <- simmap$maps[[edge_i]]
    for(segment_i in seq_along(edge_map)){
      state <- match(names(edge_map)[segment_i], state_names)
      duration <- unname(edge_map[segment_i])
      segment_a <- exp(-alpha[state] * duration)
      segment_b <- theta[state] * (1 - segment_a)
      segment_v <- sigma_sq[state] / (2 * alpha[state]) *
        (1 - exp(-2 * alpha[state] * duration))
      edge_b <- segment_a * edge_b + segment_b
      edge_v <- segment_a^2 * edge_v + segment_v
      edge_a <- segment_a * edge_a
    }
    seen <- which(assigned)
    node_mean[descendant] <- edge_a * node_mean[ancestor] + edge_b
    node_cov[descendant, seen] <- edge_a * node_cov[ancestor, seen]
    node_cov[seen, descendant] <- node_cov[descendant, seen]
    node_cov[descendant, descendant] <-
      edge_a^2 * node_cov[ancestor, ancestor] + edge_v
    assigned[descendant] <- TRUE
  }

  covariance <- node_cov[seq_len(n_tip), seq_len(n_tip), drop = FALSE]
  expected <- node_mean[seq_len(n_tip)]
  observed <- traits[match(simmap$tip.label, names(traits))]
  residual <- observed - expected
  chol_covariance <- chol(covariance)
  standardized <- backsolve(chol_covariance, residual, transpose = TRUE)
  -0.5 * (n_tip * log(2 * pi) + 2 * sum(log(diag(chol_covariance))) +
            sum(standardized^2))
}

make_random_case <- function(seed, n_tip, n_state){
  set.seed(seed)
  tree <- ape::rcoal(n_tip)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
  states <- sample(rep(seq_len(n_state), length.out = n_tip))
  names(states) <- tree$tip.label
  simmap <- phytools::make.simmap(
    tree, states, model = "ER", nsim = 1, pi = "equal", message = FALSE
  )
  traits <- stats::rnorm(n_tip, 0.5, 0.8)
  names(traits) <- tree$tip.label
  list(simmap = simmap, traits = traits)
}

make_painted_case <- function(seed, n_tip, n_state){
  set.seed(seed)
  tree <- ape::rcoal(n_tip)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
  n_all <- n_tip + tree$Nnode
  states <- integer(n_all)
  states[seq_len(n_tip)] <- sample(rep(seq_len(n_state), length.out = n_tip))
  states[(n_tip + 1L):n_all] <- sample(
    seq_len(n_state), tree$Nnode, replace = TRUE
  )
  maps <- lapply(seq_len(nrow(tree$edge)), function(i){
    from <- as.character(states[tree$edge[i, 1]])
    to <- as.character(states[tree$edge[i, 2]])
    if(from == to){
      setNames(tree$edge.length[i], from)
    }else{
      setNames(tree$edge.length[i] * c(0.4, 0.6), c(from, to))
    }
  })
  mapped_edge <- t(vapply(maps, function(edge_map){
    totals <- numeric(n_state)
    names(totals) <- as.character(seq_len(n_state))
    edge_totals <- tapply(unname(edge_map), names(edge_map), sum)
    totals[names(edge_totals)] <- edge_totals
    totals
  }, numeric(n_state)))
  colnames(mapped_edge) <- as.character(seq_len(n_state))
  tree$maps <- maps
  tree$mapped.edge <- mapped_edge
  class(tree) <- c("simmap", "phylo")
  traits <- stats::rnorm(n_tip, 0.5, 0.8)
  names(traits) <- tree$tip.label
  list(simmap = tree, traits = traits)
}

evaluate_case <- function(variant, case_name, case, parameter_set){
  simmap <- case$simmap
  traits <- case$traits
  states <- colnames(simmap$mapped.edge)
  n_state <- length(states)
  parameters <- parameter_set$values(n_state)
  alpha <- setNames(parameters$alpha, states)
  sigma_sq <- setNames(parameters$sigma_sq, states)
  theta <- setNames(parameters$theta, states)
  data <- data.frame(
    species = simmap$tip.label,
    regime = states[1],
    trait = traits[simmap$tip.label],
    stringsAsFactors = FALSE
  )

  oracle <- direct_ou_loglik(simmap, traits, alpha, sigma_sq, theta)
  likelihoods <- list()
  for(order in c("cladewise", "pruningwise")){
    ordered <- reorder_simmap(simmap, order)
    paths <- lapply(seq_len(ape::Ntip(ordered)), function(tip){
      getPathToRoot(ordered, tip)
    })
    likelihoods[[paste0("three_point_", order)]] <- OUwie.basic(
      ordered, data, simmap.tree = TRUE,
      alpha = unname(alpha), sigma.sq = unname(sigma_sq),
      theta = unname(theta), algorithm = "three.point", tip.paths = paths
    )
  }

  cladewise <- reorder_simmap(simmap, "cladewise")
  likelihoods$explicit_vcv <- OUwie.fixed(
    cladewise, data[match(cladewise$tip.label, data$species), ],
    model = "OUMVA", simmap.tree = TRUE, algorithm = "invert",
    alpha = unname(alpha), sigma.sq = unname(sigma_sq),
    theta = unname(theta), check.identify = FALSE, quiet = TRUE
  )$loglik

  data.frame(
    variant = variant, case = case_name,
    parameter_set = parameter_set$name, tips = ape::Ntip(simmap),
    states = n_state, method = names(likelihoods),
    loglik = unlist(likelihoods, use.names = FALSE),
    direct_oracle = oracle,
    abs_difference = abs(unlist(likelihoods, use.names = FALSE) - oracle),
    stringsAsFactors = FALSE
  )
}

run_worker <- function(repo, variant, output){
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
  parameter_sets <- list(
    list(name = "balanced", values = function(n) list(
      alpha = seq(0.4, 1.6, length.out = n),
      sigma_sq = seq(0.5, 1.4, length.out = n),
      theta = seq(-0.5, 2.5, length.out = n)
    )),
    list(name = "heterogeneous", values = function(n) list(
      alpha = exp(seq(log(0.08), log(5), length.out = n)),
      sigma_sq = exp(seq(log(0.15), log(3.5), length.out = n)),
      theta = seq(-2, 4, length.out = n)
    )),
    list(name = "near_bm", values = function(n) list(
      alpha = exp(seq(log(0.01), log(0.08), length.out = n)),
      sigma_sq = seq(0.2, 2, length.out = n),
      theta = seq(-1, 3, length.out = n)
    ))
  )
  case_specs <- list(
    list(name = "small_2state", seed = 1201, n_tip = 6, n_state = 2,
         generator = make_random_case),
    list(name = "medium_2state", seed = 1202, n_tip = 12, n_state = 2,
         generator = make_random_case),
    list(name = "medium_3state", seed = 1203, n_tip = 20, n_state = 3,
         generator = make_random_case),
    list(name = "larger_4state", seed = 1204, n_tip = 40, n_state = 4,
         generator = make_random_case),
    list(name = "painted_8state", seed = 1205, n_tip = 48, n_state = 8,
         generator = make_painted_case),
    list(name = "painted_16state", seed = 1206, n_tip = 64, n_state = 16,
         generator = make_painted_case)
  )
  rows <- list()
  row_i <- 0L
  for(spec in case_specs){
    case <- spec$generator(spec$seed, spec$n_tip, spec$n_state)
    for(parameters in parameter_sets){
      row_i <- row_i + 1L
      rows[[row_i]] <- evaluate_case(variant, spec$name, case, parameters)
    }
  }
  result <- do.call(rbind, rows)
  utils::write.csv(result, output, row.names = FALSE)
}

if(length(args) && args[1] == "--worker"){
  run_worker(normalizePath(args[2], mustWork = TRUE), args[3], args[4])
  quit(save = "no")
}

baseline_ref <- if(length(args) >= 1L) args[1] else "master"
current_ref <- if(length(args) >= 2L) args[2] else "HEAD"
output_dir <- file.path(dirname(script_file),
                        "houwie-upgrade-comparison-results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
snapshot_root <- tempfile("houwie-oracle-")
dir.create(snapshot_root)
on.exit(unlink(snapshot_root, recursive = TRUE), add = TRUE)
git <- Sys.which("git")
rscript <- file.path(R.home("bin"), "Rscript")
snapshot <- function(ref, name){
  archive <- file.path(snapshot_root, paste0(name, ".tar"))
  destination <- file.path(snapshot_root, name)
  dir.create(destination)
  status <- system2(git, c("-C", repo_root, "archive", ref, "-o", archive))
  if(status != 0L) stop("Could not archive Git reference ", ref)
  utils::untar(archive, exdir = destination)
  destination
}
repos <- c(master = snapshot(baseline_ref, "master"),
           current = snapshot(current_ref, "current"))
for(variant in names(repos)){
  output <- file.path(output_dir, paste0(variant, "-oracle.csv"))
  status <- system2(
    rscript, c(script_file, "--worker", repos[[variant]], variant, output)
  )
  if(status != 0L) stop("Oracle worker failed for ", variant)
}
oracle <- do.call(rbind, lapply(names(repos), function(variant){
  read.csv(file.path(output_dir, paste0(variant, "-oracle.csv")))
}))
utils::write.csv(oracle, file.path(output_dir, "oracle.csv"), row.names = FALSE)
cat("Oracle results written to ", normalizePath(output_dir), "\n", sep = "")
