#!/usr/bin/env Rscript

# A deterministic replacement for hOUwie's Monte Carlo sum over regime
# histories.
#
# hOUwie approximates
#
#   L = sum_h P(h | Q) f(x | h)
#
# by drawing a few hundred distinct histories h and adding their terms. But the
# joint model is a switching linear-Gaussian process on the tree: given the
# regime painted on an edge, the descendant's trait value is a linear-Gaussian
# function of the ancestor's. That structure means the sum can be pruned rather
# than sampled. The message a subtree sends to its parent stays a mixture of
# Gaussians in the parent's trait value; the only thing that grows is the number
# of mixture components, and that is bounded by merging components (which
# preserves total mass exactly) rather than by discarding histories (which does
# not, and is what makes the sampled sum biased).
#
# A message is
#
#   f(x) = sum_j exp(lw_j) * dnorm(x, m_j, sqrt(v_j))
#
# carried in moment form. The canonical (information) form is tempting because
# it represents a point observation as infinite precision, but evaluating a
# canonical message at a point far from its mean subtracts quantities of order
# k * x^2, and with a widely separated optimum those reach 1e19 and destroy the
# answer. Moment form never forms them: a tip observation is handled once, at
# the step where it enters an edge.

log_sum_exp <- function(x){
  finite <- x[is.finite(x)]
  if(!length(finite)) return(-Inf)
  maximum <- max(finite)
  maximum + log(sum(exp(x - maximum)))
}

# Ornstein-Uhlenbeck moments over one painted interval. The alpha -> 0 limit is
# Brownian motion and has to be taken by hand.
ou_interval <- function(alpha, sigma_sq, theta, duration){
  decay <- exp(-alpha * duration)
  variance <- if(alpha < 1e-8){
    sigma_sq * duration
  }else{
    sigma_sq * (-expm1(-2 * alpha * duration)) / (2 * alpha)
  }
  c(decay = decay, displacement = theta * (1 - decay), variance = variance)
}

# Compose two half-intervals, the first painted with `first` and the second with
# `second`, into x_descendant = A x_ancestor + B + N(0, V).
compose_half_half <- function(first, second, half, alpha, sigma_sq, theta){
  one <- ou_interval(alpha[first], sigma_sq[first], theta[first], half)
  two <- ou_interval(alpha[second], sigma_sq[second], theta[second], half)
  c(A = unname(one[["decay"]] * two[["decay"]]),
    B = unname(two[["decay"]] * one[["displacement"]] + two[["displacement"]]),
    V = unname(two[["decay"]]^2 * one[["variance"]] + two[["variance"]]))
}

# Push a message about the descendant value back through one edge, giving a
# message about the ancestor value:
#
#   integral N(x_d; A x_a + B, V) exp(g - k x_d^2 / 2 + h x_d) dx_d
#
# which is canonical in x_a with the coefficients below, where
# s = 1 / (1 + V k) is the shrinkage the edge applies. Note that k' carries a
# factor of A^2: a subtree below a branch with a strong pull says almost nothing
# about its ancestor, and this sends the message smoothly towards flat. The
# equivalent moment-form update divides by A^2 instead, so the same situation
# overflows the variance to Inf and takes the mean with it.
#
# A tip is stored as k = Inf with h holding the observed value, and enters here
# as the transition density read as a function of the ancestor.
push_through <- function(mixture, A, B, V){
  if(all(is.infinite(mixture$k))){
    value <- mixture$h
    return(list(g = mixture$g - 0.5 * log(2 * pi * V) - (value - B)^2 / (2 * V),
                h = A * (value - B) / V,
                k = rep(A^2 / V, length(mixture$g))))
  }
  g <- mixture$g; h <- mixture$h; k <- mixture$k
  s <- 1 / (1 + V * k)
  list(g = g + 0.5 * log(s) - 0.5 * k * s * B^2 + h * B * s +
         0.5 * h^2 * V * s,
       h = A * (h - k * B) * s,
       k = A^2 * k * s)
}

# Pointwise product of two mixtures: canonical parameters simply add.
multiply_mixtures <- function(left, right){
  i <- rep(seq_along(left$g), times = length(right$g))
  j <- rep(seq_along(right$g), each = length(left$g))
  list(g = left$g[i] + right$g[j],
       h = left$h[i] + right$h[j],
       k = left$k[i] + right$k[j])
}

concat_mixtures <- function(mixtures){
  list(g = unlist(lapply(mixtures, `[[`, "g"), use.names = FALSE),
       h = unlist(lapply(mixtures, `[[`, "h"), use.names = FALSE),
       k = unlist(lapply(mixtures, `[[`, "k"), use.names = FALSE))
}

empty_mixture <- function() list(g = -Inf, h = 0, k = 1)

# Reduce a mixture by repeatedly merging the pair that costs least under
# Runnalls' (2007) upper bound on the KL divergence to the unreduced mixture.
#
# Reduction is driven by `tolerance`, not by the component count. Merging stops
# as soon as the cheapest available merge would cost more than the tolerance, so
# a sharp isolated mode - which is expensive to merge by construction - survives.
# Capping the count instead is unsound: when one component dominates the
# likelihood the answer is insensitive to how many of the others are kept, so the
# value looks converged across a range of caps while a mode is being destroyed.
#
# `max_components` remains only as a hard ceiling on cost. Merges forced by it
# are performed above tolerance and their cost is recorded, so the caller can
# tell the difference between an answer that met the tolerance and one that ran
# out of budget.
#
# Merging needs moments, which exist only for a component with positive finite
# precision. Components too flat for that are pooled into one near-constant
# component rather than discarded - discarding is the mass-losing failure this
# whole approach exists to avoid.
reduce_mixture <- function(mixture, max_components = Inf, tolerance = 0,
                           diagnostics = NULL, cavity = NULL){
  n_input <- length(mixture$g)
  if(n_input <= 1L) return(mixture)
  if(!is.finite(max_components) && !(tolerance > 0)) return(mixture)

  flat <- !(mixture$k > 0) | !is.finite(mixture$k) |
    !is.finite(1 / mixture$k) | !is.finite(mixture$h)
  pooled <- NULL
  if(any(flat)){
    pooled <- list(g = log_sum_exp(mixture$g[flat]), h = 0, k = 0)
    if(!is.finite(pooled$g)) pooled <- NULL
    mixture <- list(g = mixture$g[!flat], h = mixture$h[!flat],
                    k = mixture$k[!flat])
  }
  budget <- max(max_components - as.integer(!is.null(pooled)), 1L)
  rejoin <- function(reduced){
    if(is.null(pooled)) return(reduced)
    concat_mixtures(list(reduced, pooled))
  }
  n <- length(mixture$g)
  if(!n) return(if(is.null(pooled)) empty_mixture() else pooled)
  if(n <= 1L) return(rejoin(mixture))

  variance <- 1 / mixture$k
  mean <- mixture$h / mixture$k
  log_mass <- mixture$g + 0.5 * log(2 * pi * variance) +
    0.5 * mixture$h^2 * variance
  finite_mass <- log_mass[is.finite(log_mass)]
  if(!length(finite_mass)) return(rejoin(mixture))
  offset <- max(finite_mass)
  # Moment matching has to use the component's true mass: the merged component
  # must represent the actual sum of the two it replaces.
  mass <- exp(log_mass - offset)
  mass[!is.finite(mass)] <- 0
  # Choosing *which* pair to merge is a different question, and mass is the
  # wrong currency for it. A component negligible in this message can dominate
  # once it meets its sibling at the parent, so relevance is judged against a
  # standing-in for the rest of the tree - the cavity, in expectation
  # propagation's sense - rather than against this message alone.
  weight <- if(is.null(cavity)){
    mass
  }else{
    precision <- 1 / cavity$variance + mixture$k
    shifted <- cavity$mean / cavity$variance + mixture$h
    log_relevance <- mixture$g - 0.5 * log(cavity$variance * precision) -
      cavity$mean^2 / (2 * cavity$variance) + shifted^2 / (2 * precision)
    finite_relevance <- log_relevance[is.finite(log_relevance)]
    if(!length(finite_relevance)) return(rejoin(mixture))
    scaled <- exp(log_relevance - max(finite_relevance))
    scaled[!is.finite(scaled)] <- 0
    scaled
  }

  # Runnalls' cost for merging i and j, using the identity
  #   w_i (m_i - m)^2 + w_j (m_j - m)^2 = w_i w_j (m_i - m_j)^2 / (w_i + w_j)
  # so the whole matrix can be formed without a loop.
  log_variance <- log(variance)
  pair_costs <- function(rows, cols){
    W <- outer(weight[rows], weight[cols], "+")
    merged_v <- (outer(weight[rows] * variance[rows],
                       weight[cols] * variance[cols], "+") +
                   outer(weight[rows], weight[cols], "*") *
                   outer(mean[rows], mean[cols], "-")^2 / W) / W
    0.5 * (W * log(merged_v) -
             outer(weight[rows] * log_variance[rows],
                   weight[cols] * log_variance[cols], "+"))
  }
  all_index <- seq_len(n)
  cost <- pair_costs(all_index, all_index)
  cost[!is.finite(cost)] <- Inf
  diag(cost) <- Inf
  alive <- rep(TRUE, n)
  # a per-row minimum turns the search for the cheapest pair from a scan of the
  # whole matrix into a scan of one vector, which is what makes retaining a large
  # mixture affordable at all
  row_min <- apply(cost, 1, min)
  row_arg <- max.col(-cost, ties.method = "first")
  n_alive <- n
  worst_forced <- 0

  repeat{
    if(n_alive <= 1L) break
    candidate <- which(alive)[which.min(row_min[alive])]
    best <- row_min[candidate]
    partner <- row_arg[candidate]
    if(!is.finite(best)) break
    forced <- n_alive > budget
    if(best > tolerance && !forced) break
    if(forced && best > tolerance) worst_forced <- max(worst_forced, best)
    i <- candidate
    j <- partner
    total <- mass[i] + mass[j]
    if(total <= 0) total <- .Machine$double.xmin
    centre <- (mass[i] * mean[i] + mass[j] * mean[j]) / total
    spread <- (mass[i] * variance[i] + mass[j] * variance[j] +
                 mass[i] * mass[j] * (mean[i] - mean[j])^2 / total) / total
    mass[i] <- total
    weight[i] <- weight[i] + weight[j]
    mean[i] <- centre
    variance[i] <- spread
    log_variance[i] <- log(spread)
    alive[j] <- FALSE
    n_alive <- n_alive - 1L
    cost[j, ] <- Inf
    cost[, j] <- Inf
    row_min[j] <- Inf
    others <- which(alive & seq_len(n) != i)
    if(length(others)){
      updated <- as.numeric(pair_costs(i, others))
      updated[!is.finite(updated)] <- Inf
      cost[i, others] <- updated
      cost[others, i] <- updated
    }
    # only rows whose cheapest partner just changed identity need rebuilding
    stale <- which(alive & (row_arg == i | row_arg == j))
    for(r in unique(c(i, stale))){
      if(!alive[r]) next
      row_min[r] <- min(cost[r, ])
      row_arg[r] <- which.min(cost[r, ])
    }
    if(n_alive <= 1L) break
  }

  if(!is.null(diagnostics)){
    diagnostics$forced <- max(diagnostics$forced, worst_forced)
    diagnostics$retained <- max(diagnostics$retained,
                                sum(alive) + as.integer(!is.null(pooled)))
    diagnostics$formed <- max(diagnostics$formed, n_input)
  }
  keep <- which(alive)
  precision <- 1 / variance[keep]
  rejoin(list(g = log(mass[keep]) + offset -
                0.5 * log(2 * pi * variance[keep]) -
                0.5 * mean[keep]^2 * precision,
              h = mean[keep] * precision, k = precision))
}


#' Deterministic joint likelihood for the hOUwie model.
#'
#' @param resolution regime points inserted per edge. 1 reproduces hOUwie's own
#'   midpoint history model exactly, so the result is directly comparable with
#'   an exhaustive enumeration of that model. Larger values converge on the
#'   continuous-time process at linear rather than exponential cost.
#' @param tolerance largest Runnalls merge cost allowed. Reduction stops when
#'   the cheapest available merge would cost more than this, so the approximation
#'   is controlled by an error criterion rather than by a component count. 0
#'   disables merging entirely.
#' @param max_components hard ceiling on each per-state mixture, used only to
#'   bound cost. Merges forced by it happen above `tolerance`, and the largest
#'   such cost is returned in the result's "diagnostics" attribute alongside the
#'   largest mixture retained and formed. A forced cost of 0 means every merge
#'   met the tolerance and the answer can be trusted to that level.
mixture_pruning_loglik <- function(phy, tip_states, tip_values, Q,
                                   alpha, sigma_sq, theta,
                                   root_p = NULL, resolution = 1L,
                                   tolerance = 0,
                                   max_components = 8L,
                                   edge_slack = 1,
                                   cavity_weighting = FALSE){
  phy <- ape::reorder.phylo(phy, "pruningwise")
  n_tip <- ape::Ntip(phy)
  n_states <- nrow(Q)
  if(is.null(root_p)) root_p <- rep(1 / n_states, n_states)
  tip_states <- tip_states[phy$tip.label]
  tip_values <- tip_values[phy$tip.label]

  # Reducing inside an edge happens before the message has met its sibling, so
  # it is the step with the least information about what matters. edge_slack
  # lets that reduction keep a looser bound than the node reduction, which does
  # at least see both children.
  edge_components <- if(is.finite(max_components)){
    max_components * edge_slack
  }else{
    Inf
  }
  cavity <- if(cavity_weighting){
    list(mean = mean(tip_values), variance = max(stats::var(tip_values), 1e-8))
  }else{
    NULL
  }

  diagnostics <- new.env(parent = emptyenv())
  diagnostics$forced <- 0
  diagnostics$retained <- 0L
  diagnostics$formed <- 0L

  transition_cache <- new.env(parent = emptyenv())
  transition_for <- function(duration){
    key <- sprintf("%.14g", duration)
    cached <- transition_cache[[key]]
    if(is.null(cached)){
      cached <- as.matrix(Matrix::expm(Q * duration))
      transition_cache[[key]] <- cached
    }
    cached
  }

  # Walk one edge from its descendant end to its ancestor end, one regime
  # interval at a time. `mixtures[[s]]` is always the message conditional on the
  # regime point currently being stood on being in state s.
  send_up <- function(mixtures, edge_length, allowed){
    interval <- edge_length / resolution
    transition <- transition_for(interval)
    half <- interval / 2
    composed <- lapply(seq_len(n_states), function(ancestor_state){
      lapply(seq_len(n_states), function(descendant_state){
        compose_half_half(ancestor_state, descendant_state, half,
                          alpha, sigma_sq, theta)
      })
    })
    for(step in seq_len(resolution)){
      updated <- vector("list", n_states)
      for(ancestor_state in seq_len(n_states)){
        pieces <- list()
        for(descendant_state in seq_len(n_states)){
          if(step == 1L && !is.null(allowed) &&
             !(descendant_state %in% allowed)) next
          log_step <- log(transition[ancestor_state, descendant_state])
          if(!is.finite(log_step)) next
          moments <- composed[[ancestor_state]][[descendant_state]]
          pushed <- push_through(mixtures[[descendant_state]],
                                 moments[["A"]], moments[["B"]],
                                 moments[["V"]])
          pushed$g <- pushed$g + log_step
          pieces[[length(pieces) + 1L]] <- pushed
        }
        updated[[ancestor_state]] <- if(length(pieces)){
          reduce_mixture(concat_mixtures(pieces),
                         edge_components, tolerance, diagnostics, cavity)
        }else{
          empty_mixture()
        }
      }
      mixtures <- updated
    }
    mixtures
  }

  node_mixture <- vector("list", n_tip + phy$Nnode)
  for(ancestor in unique(phy$edge[, 1])){
    combined <- NULL
    for(edge_i in which(phy$edge[, 1] == ancestor)){
      descendant <- phy$edge[edge_i, 2]
      if(descendant <= n_tip){
        # an exact observation is infinite precision, with the value parked in
        # h for push_through to read
        observation <- list(g = 0, h = unname(tip_values[descendant]),
                            k = Inf)
        child <- rep(list(observation), n_states)
        # a list entry lets one tip be compatible with several regimes, which is
        # what hidden states and ambiguous coding both need
        allowed <- if(is.list(tip_states)){
          tip_states[[descendant]]
        }else{
          unname(tip_states[descendant])
        }
      }else{
        child <- node_mixture[[descendant]]
        allowed <- NULL
      }
      message_up <- send_up(child, phy$edge.length[edge_i], allowed)
      combined <- if(is.null(combined)){
        message_up
      }else{
        lapply(seq_len(n_states), function(state){
          reduce_mixture(multiply_mixtures(combined[[state]],
                                           message_up[[state]]),
                         max_components, tolerance, diagnostics, cavity)
        })
      }
    }
    node_mixture[[ancestor]] <- combined
  }

  # hOUwie starts the continuous process at the optimum of the root's regime, so
  # the root message is read off at a point rather than integrated. Centring the
  # quadratic on that point keeps the two large terms from being formed and then
  # cancelled, which is what a widely separated optimum would otherwise do.
  root_mixture <- node_mixture[[n_tip + 1L]]
  value <- log_sum_exp(vapply(seq_len(n_states), function(state){
    mixture <- root_mixture[[state]]
    start <- theta[state]
    log(root_p[state]) +
      log_sum_exp(mixture$g + start * (mixture$h - 0.5 * mixture$k * start))
  }, numeric(1)))
  attr(value, "diagnostics") <- list(worst_forced_merge = diagnostics$forced,
                                     largest_retained = diagnostics$retained,
                                     largest_formed = diagnostics$formed)
  value
}
