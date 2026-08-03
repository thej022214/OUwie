# Deterministic joint likelihood for the hOUwie model.
#
# The sampling likelihood in hOUwie.internal.R approximates
#
#   L = sum_h P(h | Q) f(x | h)
#
# by drawing regime histories and adding their terms. The joint model is a
# switching linear-Gaussian process on the tree - given the regime painted on an
# edge, the descendant's trait value is a linear-Gaussian function of the
# ancestor's - so the same sum can be pruned instead. The message a subtree
# sends to its parent stays a mixture of Gaussians in the parent's trait value,
# and only the number of components grows.
#
# Messages are carried in canonical form
#
#   f(x) = sum_j exp(g_j - k_j x^2 / 2 + h_j x)
#
# because k = 0 is the flat function and k = Inf an exact observation, and both
# arise. Below a branch with a strong pull a subtree says almost nothing about
# its ancestor, which is k -> 0 here and a variance overflow in a mean/variance
# parameterization.
#
# The approximation is the merging of mixture components, and it is not bounded:
# see the note in reduceMessageMixture.

logSumExpSafe <- function(x){
  finite <- x[is.finite(x)]
  if(!length(finite)) return(-Inf)
  maximum <- max(finite)
  maximum + log(sum(exp(x - maximum)))
}

# Ornstein-Uhlenbeck moments over one painted interval. The alpha -> 0 limit is
# Brownian motion and has to be taken by hand.
getOUIntervalMoments <- function(alpha, sigma.sq, theta, duration){
  decay <- exp(-alpha * duration)
  variance <- if(alpha < 1e-8){
    sigma.sq * duration
  }else{
    sigma.sq * (-expm1(-2 * alpha * duration)) / (2 * alpha)
  }
  c(decay = decay, displacement = theta * (1 - decay), variance = variance)
}

# Compose two half-intervals, the first painted with `first` and the second with
# `second`, into x_descendant = A x_ancestor + B + N(0, V).
composeIntervalPair <- function(first, second, half, alpha, sigma.sq, theta){
  one <- getOUIntervalMoments(alpha[first], sigma.sq[first], theta[first], half)
  two <- getOUIntervalMoments(alpha[second], sigma.sq[second], theta[second],
                              half)
  c(A = unname(one[["decay"]] * two[["decay"]]),
    B = unname(two[["decay"]] * one[["displacement"]] + two[["displacement"]]),
    V = unname(two[["decay"]]^2 * one[["variance"]] + two[["variance"]]))
}

# Push a message about the descendant value back through one interval, giving a
# message about the ancestor value. A tip is stored as k = Inf with the observed
# value in h, and enters as the transition density read as a function of the
# ancestor.
pushMessage <- function(message, A, B, V){
  if(all(is.infinite(message$k))){
    value <- message$h
    return(list(g = message$g - 0.5 * log(2 * pi * V) - (value - B)^2 / (2 * V),
                h = A * (value - B) / V,
                k = rep(A^2 / V, length(message$g))))
  }
  shrink <- 1 / (1 + V * message$k)
  list(g = message$g + 0.5 * log(shrink) -
         0.5 * message$k * shrink * B^2 + message$h * B * shrink +
         0.5 * message$h^2 * V * shrink,
       h = A * (message$h - message$k * B) * shrink,
       k = A^2 * message$k * shrink)
}

# Pointwise product of two messages: canonical parameters simply add.
multiplyMessages <- function(left, right){
  i <- rep(seq_along(left$g), times = length(right$g))
  j <- rep(seq_along(right$g), each = length(left$g))
  list(g = left$g[i] + right$g[j],
       h = left$h[i] + right$h[j],
       k = left$k[i] + right$k[j])
}

bindMessages <- function(messages){
  list(g = unlist(lapply(messages, `[[`, "g"), use.names = FALSE),
       h = unlist(lapply(messages, `[[`, "h"), use.names = FALSE),
       k = unlist(lapply(messages, `[[`, "k"), use.names = FALSE))
}

emptyMessage <- function() list(g = -Inf, h = 0, k = 1)

# Reduce a message by repeatedly merging the pair that costs least under
# Runnalls' (2007) upper bound on the KL divergence to the unreduced mixture.
#
# This is the only approximation in the calculation and it is not bounded. On a
# tree, two sibling messages are multiplied at their parent, so a component that
# carries little weight inside one message can still dominate once it meets the
# other; the merge cost, which is local to one message, does not see that. Where
# a sharp component is merged away the likelihood can be wrong by many log units
# and, because the result is then insensitive to how many of the remaining
# components are kept, the value can look stable across a range of
# max_components while a mode is being destroyed. Treat agreement across caps as
# weak evidence, and prefer max_components = Inf whenever the tree is small
# enough to afford it.
#
# Components too flat for moments to exist are pooled into one near-constant
# component rather than discarded; discarding loses mass, which is the failure
# the sampling likelihood already has.
reduceMessageMixture <- function(message, max_components = Inf,
                                 tolerance = 0, diagnostics = NULL){
  n_input <- length(message$g)
  if(n_input <= 1L) return(message)
  if(!is.finite(max_components) && !(tolerance > 0)) return(message)

  flat <- !(message$k > 0) | !is.finite(message$k) |
    !is.finite(1 / message$k) | !is.finite(message$h)
  pooled <- NULL
  if(any(flat)){
    pooled <- list(g = logSumExpSafe(message$g[flat]), h = 0, k = 0)
    if(!is.finite(pooled$g)) pooled <- NULL
    message <- list(g = message$g[!flat], h = message$h[!flat],
                    k = message$k[!flat])
  }
  budget <- max(max_components - as.integer(!is.null(pooled)), 1L)
  rejoin <- function(reduced){
    if(is.null(pooled)) return(reduced)
    bindMessages(list(reduced, pooled))
  }
  n <- length(message$g)
  if(!n) return(if(is.null(pooled)) emptyMessage() else pooled)
  if(n <= max(budget, 1L)) return(rejoin(message))

  variance <- 1 / message$k
  centre <- message$h / message$k
  log_mass <- message$g + 0.5 * log(2 * pi * variance) +
    0.5 * message$h^2 * variance
  finite_mass <- log_mass[is.finite(log_mass)]
  if(!length(finite_mass)) return(rejoin(message))
  offset <- max(finite_mass)
  weight <- exp(log_mass - offset)
  weight[!is.finite(weight)] <- 0
  log_variance <- log(variance)

  # Runnalls' cost, using
  #   w_i (m_i - m)^2 + w_j (m_j - m)^2 = w_i w_j (m_i - m_j)^2 / (w_i + w_j)
  # so the whole matrix is formed without a loop.
  pairCosts <- function(rows, cols){
    total <- outer(weight[rows], weight[cols], "+")
    merged <- (outer(weight[rows] * variance[rows],
                     weight[cols] * variance[cols], "+") +
                 outer(weight[rows], weight[cols], "*") *
                 outer(centre[rows], centre[cols], "-")^2 / total) / total
    0.5 * (total * log(merged) -
             outer(weight[rows] * log_variance[rows],
                   weight[cols] * log_variance[cols], "+"))
  }
  index <- seq_len(n)
  cost <- pairCosts(index, index)
  cost[!is.finite(cost)] <- Inf
  diag(cost) <- Inf
  alive <- rep(TRUE, n)
  # a per-row minimum turns the search for the cheapest pair into a scan of one
  # vector, which is what makes a large retained mixture affordable
  row_min <- apply(cost, 1, min)
  row_arg <- max.col(-cost, ties.method = "first")
  n_alive <- n
  worst_forced <- 0

  while(n_alive > 1L){
    candidate <- which(alive)[which.min(row_min[alive])]
    best <- row_min[candidate]
    partner <- row_arg[candidate]
    if(!is.finite(best)) break
    forced <- n_alive > budget
    if(best > tolerance && !forced) break
    if(forced && best > tolerance) worst_forced <- max(worst_forced, best)
    combined <- weight[candidate] + weight[partner]
    if(combined <= 0) combined <- .Machine$double.xmin
    new_centre <- (weight[candidate] * centre[candidate] +
                     weight[partner] * centre[partner]) / combined
    new_variance <- (weight[candidate] * variance[candidate] +
                       weight[partner] * variance[partner] +
                       weight[candidate] * weight[partner] *
                       (centre[candidate] - centre[partner])^2 / combined) /
      combined
    weight[candidate] <- combined
    centre[candidate] <- new_centre
    variance[candidate] <- new_variance
    log_variance[candidate] <- log(new_variance)
    alive[partner] <- FALSE
    n_alive <- n_alive - 1L
    cost[partner, ] <- Inf
    cost[, partner] <- Inf
    row_min[partner] <- Inf
    others <- which(alive & index != candidate)
    if(length(others)){
      refreshed <- as.numeric(pairCosts(candidate, others))
      refreshed[!is.finite(refreshed)] <- Inf
      cost[candidate, others] <- refreshed
      cost[others, candidate] <- refreshed
    }
    stale <- which(alive & (row_arg == candidate | row_arg == partner))
    for(row_i in unique(c(candidate, stale))){
      if(!alive[row_i]) next
      row_min[row_i] <- min(cost[row_i, ])
      row_arg[row_i] <- which.min(cost[row_i, ])
    }
  }

  if(!is.null(diagnostics)){
    diagnostics$forced <- max(diagnostics$forced, worst_forced)
    diagnostics$retained <- max(diagnostics$retained,
                                sum(alive) + as.integer(!is.null(pooled)))
    diagnostics$formed <- max(diagnostics$formed, n_input)
  }
  keep <- which(alive)
  precision <- 1 / variance[keep]
  rejoin(list(g = log(weight[keep]) + offset -
                0.5 * log(2 * pi * variance[keep]) -
                0.5 * centre[keep]^2 * precision,
              h = centre[keep] * precision, k = precision))
}

# The states each tip is allowed to occupy, read off the tip rows of
# edge_liks_list so that ambiguous coding and hidden states are handled the same
# way the sampling path handles them.
getTipStateSets <- function(phy, edge_liks_list, nStates){
  n_tip <- length(phy$tip.label)
  sets <- vector("list", n_tip)
  for(edge_i in seq_len(nrow(phy$edge))){
    descendant <- phy$edge[edge_i, 2]
    if(descendant > n_tip) next
    allowed <- which(edge_liks_list[[edge_i]][1, ] > 0)
    sets[[descendant]] <- if(length(allowed)) allowed else seq_len(nStates)
  }
  sets
}

# Refuse a pruning run that cannot finish, before anything is allocated.
#
# The mixture carries one component per regime history, so its size is known in
# advance: nStates raised to the number of free regime points. On a two-state
# model that is about four million components at 24 tips and doubles with every
# further tip, and because each component is three doubles and the products at a
# node build index vectors of the full length, memory gives out well before time
# does. A run that would exceed the ceiling is stopped here rather than being
# allowed to allocate for a minute and take the session with it.
#
# The alternative - capping the mixture and merging the excess - is what this
# replaces. Its error is unbounded: a component carrying little weight in one
# message can dominate once multiplied by its sibling, so a likelihood stable
# across several ceilings is not thereby converged, and the resulting surface has
# been seen to change which model AIC selects.
checkPruningFeasible <- function(phy, nStates, resolution = 1L,
                                 max_log2_components = 26){
  regime_points <- (length(phy$tip.label) - 2L) * max(as.integer(resolution), 1L)
  log2_components <- regime_points * log2(nStates)
  if(log2_components > max_log2_components){
    stop(sprintf(paste0("algorithm = \"pruning\" is exact and its cost is fixed by the data: ",
                        "%d tips at %d states and resolution %d needs about 2^%.0f mixture components, ",
                        "past the ceiling of 2^%d this implementation will attempt. ",
                        "Use algorithm = \"sampling\", or reduce the tree to about %d tips."),
                 length(phy$tip.label), nStates, as.integer(resolution),
                 log2_components, max_log2_components,
                 floor(max_log2_components / (log2(nStates) *
                                                max(as.integer(resolution), 1L))) + 2L),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @param resolution regime points inserted per edge. 1 is hOUwie's own history
#'   model, in which the regime switches at the midpoint of an edge whose
#'   endpoints differ, and is what an exhaustive enumeration of that model
#'   reproduces. Larger values refine towards the continuous-time process at
#'   linear rather than exponential cost, first order in the step.
#' @param max_components ceiling on each per-state message. Inf is exact for the
#'   stated resolution. See reduceMessageMixture for why agreement across
#'   several finite values is not evidence of convergence.
houwiePruningLik <- function(phy, tip_state_sets, tip_values, Q,
                             alpha, sigma.sq, theta, root.p = "yang",
                             resolution = 1L, max_components = Inf,
                             tolerance = 0){
  n_tip <- length(phy$tip.label)
  nStates <- nrow(Q)
  diagnostics <- new.env(parent = emptyenv())
  diagnostics$forced <- 0
  diagnostics$retained <- 0L
  diagnostics$formed <- 0L

  transition_cache <- new.env(parent = emptyenv())
  transitionFor <- function(duration){
    key <- sprintf("%.14g", duration)
    cached <- transition_cache[[key]]
    if(is.null(cached)){
      cached <- as.matrix(expm(Q * duration))
      transition_cache[[key]] <- cached
    }
    cached
  }

  # Walk one edge from its descendant end to its ancestor end, one regime
  # interval at a time. messages[[s]] is always conditional on the regime point
  # currently stood on being in state s.
  sendUp <- function(messages, edge_length, allowed){
    interval <- edge_length / resolution
    transition <- transitionFor(interval)
    half <- interval / 2
    composed <- lapply(seq_len(nStates), function(ancestor_state){
      lapply(seq_len(nStates), function(descendant_state){
        composeIntervalPair(ancestor_state, descendant_state, half, alpha,
                            sigma.sq, theta)
      })
    })
    for(step in seq_len(resolution)){
      updated <- vector("list", nStates)
      for(ancestor_state in seq_len(nStates)){
        pieces <- list()
        for(descendant_state in seq_len(nStates)){
          if(step == 1L && !is.null(allowed) &&
             !(descendant_state %in% allowed)) next
          log_step <- log(transition[ancestor_state, descendant_state])
          if(!is.finite(log_step)) next
          moments <- composed[[ancestor_state]][[descendant_state]]
          pushed <- pushMessage(messages[[descendant_state]], moments[["A"]],
                                moments[["B"]], moments[["V"]])
          pushed$g <- pushed$g + log_step
          pieces[[length(pieces) + 1L]] <- pushed
        }
        updated[[ancestor_state]] <- if(length(pieces)){
          reduceMessageMixture(bindMessages(pieces), max_components, tolerance,
                               diagnostics)
        }else{
          emptyMessage()
        }
      }
      messages <- updated
    }
    messages
  }

  node_message <- vector("list", n_tip + phy$Nnode)
  for(ancestor in unique(phy$edge[, 1])){
    combined <- NULL
    for(edge_i in which(phy$edge[, 1] == ancestor)){
      descendant <- phy$edge[edge_i, 2]
      if(descendant <= n_tip){
        observation <- list(g = 0, h = unname(tip_values[descendant]), k = Inf)
        child <- rep(list(observation), nStates)
        allowed <- tip_state_sets[[descendant]]
      }else{
        child <- node_message[[descendant]]
        allowed <- NULL
      }
      message_up <- sendUp(child, phy$edge.length[edge_i], allowed)
      combined <- if(is.null(combined)){
        message_up
      }else{
        lapply(seq_len(nStates), function(state){
          reduceMessageMixture(multiplyMessages(combined[[state]],
                                                message_up[[state]]),
                               max_components, tolerance, diagnostics)
        })
      }
    }
    node_message[[ancestor]] <- combined
  }

  # hOUwie starts the continuous process at the optimum of the root's regime, so
  # the root message is read off at a point. Centring the quadratic on that point
  # avoids forming and then cancelling two large terms.
  root_message <- node_message[[n_tip + 1L]]
  per_state <- vapply(seq_len(nStates), function(state){
    message <- root_message[[state]]
    start <- theta[state]
    logSumExpSafe(message$g + start * (message$h - 0.5 * message$k * start))
  }, numeric(1))

  root_liks <- if(inherits(root.p[1], what = "character")){
    if(root.p == "yang"){
      stationary <- c(MASS::Null(Q))
      stationary / sum(stationary)
    }else if(root.p == "maddfitz"){
      weights <- exp(per_state - max(per_state[is.finite(per_state)]))
      weights / sum(weights)
    }else{
      rep(1 / nStates, nStates)
    }
  }else{
    root.p / sum(root.p)
  }
  value <- logSumExpSafe(log(root_liks) + per_state)
  attr(value, "diagnostics") <- list(worst_forced_merge = diagnostics$forced,
                                     largest_retained = diagnostics$retained,
                                     largest_formed = diagnostics$formed)
  value
}
