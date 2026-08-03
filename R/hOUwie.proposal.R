# Continuous-trait factors for the history proposal.
#
# The optimal proposal for a regime history is q*(h) = P(h | Q, discrete, x), and
# sequential importance sampling says how to approach it: weight each state at
# each node by the partial likelihood of all the data below it. For the discrete
# data that partial likelihood is what getConditionalInternodeLik already
# propagates. The continuous analogue is not a number but a function of the trait
# value y at that node,
#
#   D_s(y) = p(continuous data below | regime = s here, trait = y here)
#
# which is the message a pruning pass carries. Collapsed to one Gaussian per node
# per state it is cheap and linear in the tree, and because this only ever builds
# a proposal, the collapse costs effective sample size and never biases the
# likelihood.
#
# Two things distinguish this from getCherryConditionals, which it replaces.
# Descendants combine by multiplying messages rather than by averaging
# normalized per-tip vectors, so evidence accumulates instead of washing out. And
# each tip's evidence enters exactly once, propagated through the real tree,
# rather than being injected at every ancestor on its path to the root and then
# propagated again - two errors that partly cancelled and made the old version
# impossible to reason about.

# Integral of a canonical-form message against a reference density for the trait
# at the node, which turns the message into the scalar the sampler needs. The
# reference is the stationary distribution of the regime being scored: it asks
# how well the data below fit if this node sits in state s. A message from a tip
# is a point observation and integrates to the density at that point.
logMessageAgainstReference <- function(message, mean_ref, var_ref){
  if(!is.finite(message$k)){
    return(message$g + stats::dnorm(message$h, mean_ref, sqrt(var_ref),
                                    log = TRUE))
  }
  precision <- 1 / var_ref + message$k
  centre <- (mean_ref / var_ref + message$h) / precision
  message$g - 0.5 * log(1 + message$k * var_ref) +
    0.5 * precision * centre^2 - mean_ref^2 / (2 * var_ref)
}

# Stationary variance of each regime. Under weak pull the stationary
# distribution does not exist, so the reference falls back to something diffuse
# on the scale of the data; any bounded choice is legitimate here because the
# reference only shapes the proposal.
getReferenceVariances <- function(alpha, sigma.sq, diffuse){
  vapply(seq_along(alpha), function(state){
    if(alpha[state] < 1e-8) return(diffuse)
    min(sigma.sq[state] / (2 * alpha[state]), diffuse)
  }, numeric(1))
}

#' @param temper exponent applied to the factors. The reference is the
#'   stationary distribution of a regime, which asserts that a lineage has
#'   already equilibrated; a lineage that switched recently still carries its old
#'   trait value, and under a sharp reference such a history is proposed almost
#'   never and carries an enormous weight when it is. Flattening the factor is
#'   the standard remedy - it costs discrimination and buys bounded weights, and
#'   because this only shapes the proposal it cannot bias anything either way.
#' @return a matrix of log factors, one row per node (tips first, as in the
#'   phylo node numbering) and one column per state.
getContinuousNodeFactors <- function(phy, tip_values, Q, alpha, sigma.sq, theta,
                                     temper = 1){
  n_tip <- length(phy$tip.label)
  n_state <- nrow(Q)
  diffuse <- max(stats::var(tip_values, na.rm = TRUE) * 4, 1e-6)
  var_ref <- getReferenceVariances(alpha, sigma.sq, diffuse)

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

  node_message <- vector("list", n_tip + phy$Nnode)
  log_factors <- matrix(0, n_tip + phy$Nnode, n_state)

  for(tip in seq_len(n_tip)){
    observation <- list(g = 0, h = unname(tip_values[tip]), k = Inf)
    node_message[[tip]] <- rep(list(observation), n_state)
    log_factors[tip, ] <- vapply(seq_len(n_state), function(state)
      logMessageAgainstReference(observation, theta[state], var_ref[state]),
      numeric(1))
  }

  # phy arrives in pruningwise order, so ancestors are met only once every
  # descendant of theirs has a message
  for(ancestor in unique(phy$edge[, 1])){
    combined <- NULL
    for(edge_i in which(phy$edge[, 1] == ancestor)){
      descendant <- phy$edge[edge_i, 2]
      duration <- phy$edge.length[edge_i]
      transition <- transitionFor(duration)
      half <- duration / 2
      from_child <- lapply(seq_len(n_state), function(ancestor_state){
        pieces <- lapply(seq_len(n_state), function(descendant_state){
          moments <- composeIntervalPair(ancestor_state, descendant_state, half,
                                         alpha, sigma.sq, theta)
          pushed <- pushMessage(node_message[[descendant]][[descendant_state]],
                                moments[["A"]], moments[["B"]], moments[["V"]])
          pushed$g <- pushed$g + log(transition[ancestor_state,
                                                descendant_state])
          pushed
        })
        # one Gaussian per state is the whole point: the mixture is what makes an
        # exact pruning likelihood unaffordable, and a proposal does not need it
        reduceMessageMixture(bindMessages(pieces), max_components = 1L)
      })
      combined <- if(is.null(combined)){
        from_child
      }else{
        lapply(seq_len(n_state), function(state)
          reduceMessageMixture(multiplyMessages(combined[[state]],
                                                from_child[[state]]),
                               max_components = 1L))
      }
    }
    node_message[[ancestor]] <- combined
    log_factors[ancestor, ] <- vapply(seq_len(n_state), function(state)
      logMessageAgainstReference(combined[[state]], theta[state],
                                 var_ref[state]), numeric(1))
  }
  log_factors * temper
}

# Multiply the factors into the already-propagated discrete conditionals. This
# happens after getConditionalInternodeLik rather than before it, which is what
# keeps each tip's continuous evidence from being counted once per ancestor: the
# messages have already done their own propagation, so letting the discrete pass
# propagate them again would count the same data twice.
applyContinuousNodeFactors <- function(conditional_probs, phy, log_factors){
  edge_liks_list <- conditional_probs$edge_liks_list
  n_tip <- length(phy$tip.label)
  scaleRow <- function(row, node){
    weights <- row * exp(log_factors[node, ] - max(log_factors[node, ]))
    total <- sum(weights)
    if(!is.finite(total) || total <= 0) return(row)
    weights / total
  }
  for(edge_i in seq_len(nrow(phy$edge))){
    liks <- edge_liks_list[[edge_i]]
    # row 1 is the descendant end of the edge and the last row the ancestor end;
    # at one regime point per edge those are exactly the two tree nodes
    liks[1, ] <- scaleRow(liks[1, ], phy$edge[edge_i, 2])
    liks[nrow(liks), ] <- scaleRow(liks[nrow(liks), ], phy$edge[edge_i, 1])
    edge_liks_list[[edge_i]] <- liks
  }
  conditional_probs$edge_liks_list <- edge_liks_list
  root <- n_tip + 1L
  conditional_probs$root_state <- scaleRow(conditional_probs$root_state, root)
  conditional_probs
}
