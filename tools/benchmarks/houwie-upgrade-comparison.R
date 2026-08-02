#!/usr/bin/env Rscript

# Reproducible comparison of the houwie-update branch with the repository's
# primary branch. Each revision is loaded in a separate R process so package
# namespaces and compiled functions cannot leak between variants.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."),
                           mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)

map_signature <- function(fit){
  lapply(fit$simmaps, function(tree){
    lapply(tree$maps, function(edge){
      list(state = names(edge), duration = unname(edge))
    })
  })
}

map_fingerprint <- function(simmap){
  edge_order <- order(simmap$edge[, 1], simmap$edge[, 2])
  edge_text <- vapply(edge_order, function(i){
    segments <- paste(
      names(simmap$maps[[i]]),
      formatC(unname(simmap$maps[[i]]), digits = 17, format = "fg"),
      sep = ":", collapse = ","
    )
    paste(simmap$edge[i, 1], simmap$edge[i, 2], segments, sep = ">")
  }, character(1))
  paste(edge_text, collapse = "|")
}

# Independent linear-Gaussian calculation for a piecewise OU process. It does
# not use OUwie's weight matrix, transformed tree, covariance, or three-point
# implementation. hOUwie fixes the root at its regime optimum by default, so
# the root contributes no variance here.
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
  assigned[root] <- TRUE

  edge_order <- ape::reorder.phylo(simmap, "cladewise", index.only = TRUE)
  for(edge_i in edge_order){
    ancestor <- simmap$edge[edge_i, 1]
    descendant <- simmap$edge[edge_i, 2]
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

log_sum_exp <- function(x){
  max(x) + log(sum(exp(x - max(x))))
}

make_case <- function(label, n_tips, rate_cat, n_sim, seed, repetitions){
  set.seed(seed)
  phy <- ape::rcoal(n_tips)
  data <- data.frame(
    species = phy$tip.label,
    state = rep(1:2, length.out = n_tips),
    trait = 1.5 + seq(-0.4, 0.4, length.out = n_tips) +
      stats::rnorm(n_tips, 0, 0.08),
    stringsAsFactors = FALSE
  )
  hdat <- organizeHOUwieDat(data, "none", TRUE)
  index_disc <- getDiscreteModel(data[, 1:2], "ER", rate_cat, FALSE, TRUE)
  index_disc[index_disc == 0] <- NA
  # Use only helpers present on both revisions. The benchmark data are encoded
  # as consecutive integer states, so the maximum is the observed state count.
  n_states <- max(as.numeric(hdat$data.cor[, 2]), na.rm = TRUE)
  index_cont <- getOUParamStructure("OUM", n_states, rate_cat, TRUE)
  counts <- vapply(seq_len(3), function(i){
    length(unique(stats::na.omit(index_cont[i, ])))
  }, integer(1))
  p <- c(
    rep(0.2, max(index_disc, na.rm = TRUE)),
    rep(0.5, counts[1]),
    rep(0.8, counts[2]),
    seq(1.25, 1.75, length.out = counts[3])
  )
  list(label = label, phy = phy, data = data, p = p, nSim = n_sim,
       rate.cat = rate_cat, seed = seed, states = nrow(index_disc),
       parameters = length(p), repetitions = repetitions)
}

fit_case <- function(case, seed, legacy_rng = FALSE){
  fit_args <- list(
    phy = case$phy,
    data = case$data,
    rate.cat = case$rate.cat,
    discrete_model = "ER",
    continuous_model = "OUM",
    null.model = TRUE,
    nSim = case$nSim,
    p = case$p,
    root.p = "flat",
    quiet = TRUE,
    sample_tips = FALSE,
    sample_nodes = FALSE,
    adaptive_sampling = FALSE
  )
  if(legacy_rng && "common_random_numbers" %in% names(formals(hOUwie))){
    fit_args$common_random_numbers <- FALSE
  }
  set.seed(seed)
  suppressWarnings(do.call(hOUwie, fit_args))
}

safe_check <- function(name, code, threshold = NA_real_, detail = ""){
  tryCatch({
    result <- force(code)
    data.frame(check = name, passed = isTRUE(result$passed),
               value = as.numeric(result$value), threshold = threshold,
               detail = detail, stringsAsFactors = FALSE)
  }, error = function(e){
    data.frame(check = name, passed = FALSE, value = NA_real_,
               threshold = threshold, detail = conditionMessage(e),
               stringsAsFactors = FALSE)
  })
}

run_checks <- function(){
  data(tworegime, envir = environment())

  make_simmap <- function(){
    map <- getMapFromNode(tree, trait[, 2], tree$node.label, 0.5)
    getMapFromSubstHistory(list(map), tree)[[1]]
  }
  reorder_simmap <- function(simmap, order){
    index <- reorder.phylo(simmap, order, index.only = TRUE)
    simmap$edge <- simmap$edge[index, ]
    simmap$edge.length <- simmap$edge.length[index]
    simmap$maps <- simmap$maps[index]
    simmap$mapped.edge <- simmap$mapped.edge[index, ]
    attr(simmap, "order") <- order
    simmap
  }

  checks <- list()
  checks[[length(checks) + 1L]] <- safe_check(
    "edge-order invariant continuous likelihood",
    {
      simmap <- make_simmap()
      data.ou <- data.frame(sp = trait[, 1], reg = trait[, 2], x = trait[, 3])
      p <- list(alpha = c(2, 4), sigma.sq = c(1, 0.5), theta = c(1, 2))
      likelihoods <- vapply(c("cladewise", "pruningwise"), function(order){
        phy <- reorder_simmap(simmap, order)
        paths <- lapply(seq_len(Ntip(phy)), function(tip) getPathToRoot(phy, tip))
        OUwie.basic(phy, data.ou, simmap.tree = TRUE, alpha = p$alpha,
                    sigma.sq = p$sigma.sq, theta = p$theta,
                    algorithm = "three.point", tip.paths = paths)
      }, numeric(1))
      delta <- abs(diff(likelihoods))
      list(passed = delta <= 1e-8, value = delta)
    }, threshold = 1e-8,
    detail = "Same painted tree in cladewise and pruningwise edge order"
  )

  checks[[length(checks) + 1L]] <- safe_check(
    "mapped-edge state-column order is invariant",
    {
      simmap <- make_simmap()
      data.ou <- data.frame(sp = trait[, 1], reg = trait[, 2], x = trait[, 3])
      evaluate <- function(current.map){
        paths <- lapply(seq_len(Ntip(current.map)), function(tip){
          getPathToRoot(current.map, tip)
        })
        OUwie.basic(
          current.map, data.ou, simmap.tree = TRUE,
          alpha = c(0.5, 0.5), sigma.sq = c(0.8, 0.8),
          theta = c(1, 3), algorithm = "three.point", tip.paths = paths
        )
      }
      original <- evaluate(simmap)
      permuted <- simmap
      permuted$mapped.edge <- permuted$mapped.edge[
        , rev(seq_len(ncol(permuted$mapped.edge))), drop = FALSE
      ]
      delta <- abs(original - evaluate(permuted))
      list(passed = delta <= 1e-8, value = delta)
    }, threshold = 1e-8,
    detail = paste(
      "Permuting mapped.edge columns must not reassign alpha, sigma2, or theta",
      "to different biological states"
    )
  )

  checks[[length(checks) + 1L]] <- safe_check(
    "three-point likelihood agrees with explicit VCV",
    {
      simmap <- make_simmap()
      data.ou <- data.frame(sp = trait[, 1], reg = trait[, 2], x = trait[, 3])
      alpha <- c(2, 4); sigma.sq <- c(1, 0.5); theta <- c(1, 2)
      pruning <- reorder_simmap(simmap, "pruningwise")
      paths <- lapply(seq_len(Ntip(pruning)), function(tip)
        getPathToRoot(pruning, tip))
      three_point <- OUwie.basic(
        pruning, data.ou, simmap.tree = TRUE, alpha = alpha,
        sigma.sq = sigma.sq, theta = theta, algorithm = "three.point",
        tip.paths = paths
      )
      cladewise <- reorder_simmap(simmap, "cladewise")
      explicit <- OUwie.fixed(
        cladewise,
        data.ou[match(cladewise$tip.label, data.ou$sp), ],
        model = "OUMVA", simmap.tree = TRUE, algorithm = "invert",
        alpha = alpha, sigma.sq = sigma.sq, theta = theta,
        check.identify = FALSE, quiet = TRUE
      )$loglik
      delta <- abs(as.numeric(three_point) - as.numeric(explicit))
      list(passed = delta <= 1e-6, value = delta)
    }, threshold = 1e-6,
    detail = "Independent covariance-matrix likelihood is the reference"
  )

  checks[[length(checks) + 1L]] <- safe_check(
    "likelihood cache uses the complete parameter vector",
    {
      simmap <- make_simmap()
      hdat <- organizeHOUwieDat(trait, "none", TRUE)
      index.disc <- getDiscreteModel(hdat$data.cor, "ER", 1, FALSE, TRUE)
      index.disc[index.disc == 0] <- NA
      index.cont <- getOUParamStructure("OUMV", 2, 1, FALSE)
      paths <- lapply(seq_len(Nnode(simmap) + Ntip(simmap)), function(node)
        getPathToRoot(simmap, node))
      edge.liks <- getEdgeLiks(simmap, hdat$data.cor, 2, 1,
                               max(branching.times(simmap)) + 1)
      n.p <- max(index.disc, na.rm = TRUE) + max(index.cont, na.rm = TRUE)
      empty_cache <- function(){
        data.table::as.data.table(data.frame(matrix(
          c(0, rep(1e5, n.p)), byrow = TRUE, ncol = n.p + 1, nrow = 50
        )))
      }
      likelihood <- function(p, cache){
        hOUwie.fixed.dev(
          p = log(p), simmaps = list(simmap), data = hdat$data.ou,
          rate.cat = 1, tip.fog = "none", index.disc = index.disc,
          index.cont = index.cont, root.p = "yang",
          edge_liks_list = edge.liks, all.paths = paths,
          split.liks = FALSE, global_liks_mat = cache
        )
      }
      pars <- c(1, 2, 0.5, 4, 1, 3)
      swapped <- c(1, 2, 4, 0.5, 1, 3)
      shared <- empty_cache()
      likelihood(pars, shared)
      cached <- likelihood(swapped, shared)
      fresh <- likelihood(swapped, empty_cache())
      delta <- abs(cached - fresh)
      list(passed = delta <= 1e-8, value = delta)
    }, threshold = 1e-8,
    detail = "Vectors with the same sum but different sigma values must not collide"
  )

  checks[[length(checks) + 1L]] <- safe_check(
    "tip sampling is invariant to data-row order",
    {
      set.seed(42)
      phy <- rcoal(8)
      phy$edge.length <- phy$edge.length / max(branching.times(phy))
      dat <- data.frame(
        sp = phy$tip.label,
        reg = c(1, 1, 2, 2, 1, 2, 1, 2),
        x = c(1.2, 0.9, 3.1, 3.4, 1.0, 2.9, 1.1, 3.2)
      )
      p <- c(0.5, 0.5, 0.3, 0.3, 1, 1, 1.5, 3, 2, 2.5)
      evaluate <- function(current.data){
        fit.args <- list(
          phy = phy, data = current.data, rate.cat = 2,
          discrete_model = "ER", continuous_model = "OUM", nSim = 10,
          p = p, sample_tips = TRUE, quiet = TRUE
        )
        if("common_random_numbers" %in% names(formals(hOUwie))){
          fit.args$common_random_numbers <- FALSE
        }
        set.seed(7)
        do.call(hOUwie, fit.args)$loglik
      }
      set.seed(1)
      shuffled <- dat[sample(nrow(dat)), ]
      delta <- abs(evaluate(shuffled) - evaluate(dat))
      list(passed = delta <= 1e-8, value = delta)
    }, threshold = 1e-8,
    detail = "Species names, not row or pruningwise edge positions, pair observations"
  )

  checks[[length(checks) + 1L]] <- safe_check(
    "common-random-number objective is available",
    {
      available <- exists("makeCommonRandomObjective", mode = "function",
                          inherits = TRUE) &&
        "common_random_numbers" %in% names(formals(hOUwie))
      if(!available){
        list(passed = FALSE, value = 0)
      }else{
        set.seed(99)
        caller.seed <- .Random.seed
        objective <- makeCommonRandomObjective(function(p) p + runif(4), 123)
        first <- objective(1)
        second <- objective(2)
        stable <- identical(.Random.seed, caller.seed) &&
          isTRUE(all.equal(second - first, rep(1, 4)))
        list(passed = stable, value = as.numeric(stable))
      }
    }, threshold = 1,
    detail = "Repeated comparisons reuse draws without leaking RNG state"
  )

  checks[[length(checks) + 1L]] <- safe_check(
    "multi-digit state histories have collision-free identifiers",
    {
      available <- exists("stateSampleId", mode = "function", inherits = TRUE)
      distinct <- available &&
        !identical(stateSampleId(list(c(1, 11))),
                   stateSampleId(list(c(11, 1))))
      list(passed = distinct, value = as.numeric(distinct))
    }, threshold = 1,
    detail = "Histories (1,11) and (11,1) must remain distinct"
  )

  do.call(rbind, checks)
}

if(length(args) && args[1] == "--worker"){
  repo <- normalizePath(args[2], mustWork = TRUE)
  variant <- args[3]
  output_dir <- args[4]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
  data.table::setDTthreads(1L)

  cases <- list(
    make_case("moderate_nSim25", 96L, 2L, 25L, 8101L, 5L),
    make_case("moderate_nSim100", 96L, 2L, 100L, 8102L, 5L),
    make_case("tree_1000", 1000L, 2L, 25L, 8103L, 3L),
    make_case("states_16", 250L, 8L, 25L, 8104L, 3L),
    make_case("combined", 500L, 4L, 50L, 8105L, 3L)
  )
  invisible(fit_case(make_case("warmup", 30L, 2L, 5L, 8001L, 1L),
                     8001L))

  benchmark <- do.call(rbind, lapply(cases, function(case){
    do.call(rbind, lapply(seq_len(case$repetitions), function(repetition){
      gc(FALSE)
      elapsed <- system.time({
        fit <- fit_case(case, case$seed + repetition)
      })[["elapsed"]]
      data.frame(
        variant = variant, case = case$label, repetition = repetition,
        tips = Ntip(case$phy), nSim = case$nSim, states = case$states,
        parameters = case$parameters, elapsed_seconds = unname(elapsed),
        loglik = fit$loglik, maps_returned = length(fit$simmaps),
        map_bytes = as.numeric(object.size(fit$simmaps)),
        stringsAsFactors = FALSE
      )
    }))
  }))
  utils::write.csv(benchmark,
                   file.path(output_dir, paste0(variant, "-benchmark.csv")),
                   row.names = FALSE)

  checks <- run_checks()
  checks$variant <- variant
  checks <- checks[, c("variant", setdiff(names(checks), "variant"))]
  utils::write.csv(checks,
                   file.path(output_dir, paste0(variant, "-checks.csv")),
                   row.names = FALSE)

  equivalence.case <- make_case("equivalence", 96L, 2L, 25L, 8111L, 1L)
  equivalence.fit <- fit_case(equivalence.case, equivalence.case$seed,
                              legacy_rng = TRUE)
  state_names <- as.character(seq_len(ncol(equivalence.fit$solution.cont)))
  alpha <- setNames(unname(equivalence.fit$solution.cont[1, ]), state_names)
  sigma_sq <- setNames(unname(equivalence.fit$solution.cont[2, ]), state_names)
  theta <- setNames(unname(equivalence.fit$solution.cont[3, ]), state_names)
  traits <- setNames(equivalence.case$data$trait,
                     equivalence.case$data$species)
  independent_continuous <- vapply(equivalence.fit$simmaps, function(simmap){
    direct_ou_loglik(simmap, traits, alpha, sigma_sq, theta)
  }, numeric(1))
  map_rows <- data.frame(
    fingerprint = vapply(equivalence.fit$simmaps, map_fingerprint,
                         character(1)),
    discrete = equivalence.fit$all_disc_liks,
    stored_continuous = equivalence.fit$all_cont_liks,
    independent_continuous = independent_continuous,
    stringsAsFactors = FALSE
  )
  saveRDS(list(
    loglik = equivalence.fit$loglik,
    discrete = equivalence.fit$all_disc_liks,
    continuous = equivalence.fit$all_cont_liks,
    maps = map_signature(equivalence.fit),
    map_rows = map_rows,
    independent_total = log_sum_exp(map_rows$discrete +
                                      map_rows$independent_continuous)
  ), file.path(output_dir, paste0(variant, "-equivalence.rds")))
  quit(save = "no")
}

if(length(args) && args[1] == "--tests"){
  repo <- normalizePath(args[2], mustWork = TRUE)
  variant <- args[3]
  output <- args[4]
  result <- try(devtools::test(repo, reporter = "summary"), silent = TRUE)
  if(inherits(result, "try-error")){
    passed <- FALSE
    detail <- as.character(result)
  } else {
    summary <- as.data.frame(result)
    passed <- all(summary$failed == 0L & !summary$error)
    detail <- if(passed){
      paste(sum(summary$passed), "expectations passed across", nrow(summary), "tests")
    } else {
      paste(sum(summary$failed), "failed expectations and", sum(summary$error), "errors")
    }
  }
  utils::write.csv(data.frame(
    variant = variant, passed = passed,
    detail = detail,
    stringsAsFactors = FALSE
  ), output, row.names = FALSE)
  quit(save = "no", status = if(passed) 0L else 1L)
}

baseline_ref <- if(length(args) >= 1L) args[1] else "master"
current_ref <- if(length(args) >= 2L) args[2] else "HEAD"
output_dir <- file.path(dirname(script_file),
                        "houwie-upgrade-comparison-results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
snapshot_root <- tempfile("houwie-branch-comparison-")
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
refs <- data.frame(
  variant = names(repos),
  requested_ref = c(baseline_ref, current_ref),
  commit = vapply(c(baseline_ref, current_ref), function(ref){
    system2(git, c("-C", repo_root, "rev-parse", ref), stdout = TRUE)
  }, character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(refs, file.path(output_dir, "refs.csv"), row.names = FALSE)

for(variant in names(repos)){
  status <- system2(rscript,
                    c(script_file, "--worker", repos[[variant]], variant,
                      output_dir))
  if(status != 0L) stop("Benchmark worker failed for ", variant)
}

benchmark <- do.call(rbind, lapply(names(repos), function(variant){
  read.csv(file.path(output_dir, paste0(variant, "-benchmark.csv")))
}))
checks <- do.call(rbind, lapply(names(repos), function(variant){
  read.csv(file.path(output_dir, paste0(variant, "-checks.csv")))
}))
utils::write.csv(benchmark, file.path(output_dir, "benchmark.csv"),
                 row.names = FALSE)
utils::write.csv(checks, file.path(output_dir, "checks.csv"), row.names = FALSE)

baseline <- readRDS(file.path(output_dir, "master-equivalence.rds"))
current <- readRDS(file.path(output_dir, "current-equivalence.rds"))
max_difference <- function(x, y){
  if(length(x) != length(y)) return(Inf)
  max(abs(as.numeric(x) - as.numeric(y)))
}
matched_maps <- merge(
  baseline$map_rows, current$map_rows, by = "fingerprint",
  suffixes = c("_master", "_current")
)
equivalence <- data.frame(
  maps_identical = identical(baseline$maps, current$maps),
  same_history_set = setequal(baseline$map_rows$fingerprint,
                              current$map_rows$fingerprint),
  histories_matched = nrow(matched_maps),
  discrete_max_abs_difference = max_difference(
    matched_maps$discrete_master, matched_maps$discrete_current
  ),
  independent_continuous_max_abs_difference = max_difference(
    matched_maps$independent_continuous_master,
    matched_maps$independent_continuous_current
  ),
  master_stored_continuous_max_abs_error = max(abs(
    baseline$map_rows$stored_continuous -
      baseline$map_rows$independent_continuous
  )),
  current_stored_continuous_max_abs_error = max(abs(
    current$map_rows$stored_continuous -
      current$map_rows$independent_continuous
  )),
  master_reported_total = baseline$loglik,
  current_reported_total = current$loglik,
  master_independent_total = baseline$independent_total,
  current_independent_total = current$independent_total,
  reported_total_abs_difference = abs(baseline$loglik - current$loglik),
  master_reported_abs_error = abs(baseline$loglik -
                                    baseline$independent_total),
  current_reported_abs_error = abs(current$loglik -
                                     current$independent_total)
)
utils::write.csv(equivalence, file.path(output_dir, "equivalence.csv"),
                 row.names = FALSE)

test_rows <- list()
for(variant in names(repos)){
  test_output <- file.path(output_dir, paste0(variant, "-tests.csv"))
  log_output <- file.path(output_dir, paste0(variant, "-tests.log"))
  status <- system2(rscript,
                    c(script_file, "--tests", repos[[variant]], variant,
                      test_output), stdout = log_output, stderr = log_output)
  if(!file.exists(test_output)){
    utils::write.csv(data.frame(
      variant = variant, passed = FALSE,
      detail = paste("Test worker exited with status", status)
    ), test_output, row.names = FALSE)
  }
  test_rows[[variant]] <- read.csv(test_output)
}
utils::write.csv(do.call(rbind, test_rows),
                 file.path(output_dir, "tests.csv"), row.names = FALSE)

cat("Comparison results written to ", normalizePath(output_dir), "\n", sep = "")
