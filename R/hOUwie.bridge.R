# Exact endpoint-conditioned regime paths by uniformization.
#
# hOUwie's own history model paints an edge with two halves, so a regime shift
# can only ever sit at the midpoint of a branch. That is an approximation to the
# continuous-time process, and refining it by slicing edges more finely is a bad
# trade: every extra slice multiplies the size of the discrete history space the
# proposal has to cover, and the effective sample size falls faster than the
# painting improves.
#
# Drawing the within-branch path from the true CTMC bridge avoids that. Sampling
# a history in two stages - node states from the discrete posterior, then the
# path along each edge conditional on its two endpoints - makes the proposal
#
#   q(h) = q(node states) * prod_e [ pathdensity(h_e | Q) / P_ab(t_e) ]
#
# while the model density is root * prod_e pathdensity(h_e | Q). The path
# densities cancel, so the importance weight is still the node-state quantity
# getStateSampleProb already computes and nothing about the discrete calculation
# changes. The cancellation is what makes this exact, and it holds only if the
# bridge itself is an exact draw - which is why rejection sampling with a
# fallback path is not usable here.
#
# Uniformization (Hobolth & Stone 2009): with mu = max_i(-Q_ii) and
# R = I + Q / mu, a CTMC path is a Poisson(mu * t) number of virtual events whose
# states follow the discrete chain R, some of which are self-transitions.
# Conditioning on the endpoints keeps that structure, so the bridge can be drawn
# without rejection and always terminates. R and its powers do not depend on the
# branch, so they are built once per parameter evaluation and every branch after
# that is a table lookup and some uniform draws; that amortization is what
# rejection sampling cannot do.

# Powers of R up to the largest event count any branch will plausibly need. NULL
# means the count needed is past the cap, which only happens where the rates are
# large enough that the likelihood is not meaningful anyway; the caller treats
# that as a failed evaluation rather than truncating the draw silently.
makeUniformizationTable <- function(Q, max_branch_length, max_events_cap = 2000L){
  n_state <- nrow(Q)
  mu <- max(-diag(Q))
  if(!is.finite(mu) || mu <= 0){
    return(list(mu = 0, R = diag(n_state),
                powers = array(diag(n_state), c(n_state, n_state, 1L)),
                max_events = 0L))
  }
  n_needed <- as.integer(stats::qpois(1 - 1e-12, mu * max_branch_length)) + 10L
  if(!is.finite(n_needed) || n_needed > max_events_cap) return(NULL)
  R <- diag(n_state) + Q / mu
  powers <- array(0, c(n_state, n_state, n_needed + 1L))
  powers[, , 1L] <- diag(n_state)
  for(i in seq_len(n_needed)){
    powers[, , i + 1L] <- powers[, , i] %*% R
  }
  list(mu = mu, R = R, powers = powers, max_events = n_needed)
}

# One branch, conditional on the states at both ends. Returns segment durations
# named by state, ordered rootward to tipward, or NULL if the draw is not
# representable.
sampleBranchBridge <- function(init, final, blen, uni){
  if(blen <= 0){
    return(stats::setNames(blen, init))
  }
  if(uni$mu <= 0){
    if(init != final) return(NULL)
    return(stats::setNames(blen, init))
  }
  n_grid <- 0:uni$max_events
  # P(n events | endpoints) is the Poisson weight times R^n[a, b]; normalizing
  # over the grid avoids having to divide by P_ab(t) and keeps this consistent
  # with itself rather than with a separately computed matrix exponential
  weights <- stats::dpois(n_grid, uni$mu * blen) *
    uni$powers[init, final, n_grid + 1L]
  weights[!is.finite(weights)] <- 0
  if(sum(weights) <= 0) return(NULL)
  n_events <- sample(n_grid, 1L, prob = weights)
  if(n_events == 0L){
    # no events at all is only consistent with the endpoints agreeing, which the
    # weights above already enforce because R^0 is the identity
    return(stats::setNames(blen, init))
  }
  times <- sort(stats::runif(n_events, 0, blen))
  states <- integer(n_events + 1L)
  states[1L] <- init
  states[n_events + 1L] <- final
  if(n_events > 1L){
    for(i in seq_len(n_events - 1L)){
      # the state after event i is drawn from one step of R weighted by the
      # chance of still reaching `final` in the events that remain
      step <- uni$R[states[i], ] * uni$powers[, final, n_events - i + 1L]
      step[!is.finite(step)] <- 0
      if(sum(step) <= 0) return(NULL)
      states[i + 1L] <- sample.int(length(step), 1L, prob = step)
    }
  }
  durations <- diff(c(0, times, blen))
  # self-transitions are the virtual events uniformization introduces and are not
  # part of the path, so runs of one state collapse to a single segment
  runs <- rle(states)
  grouping <- rep(seq_along(runs$lengths), runs$lengths)
  stats::setNames(as.numeric(tapply(durations, grouping, sum)), runs$values)
}

# Paint a whole sampled history. state_sample[[edge]] holds the states at the
# regime points of that edge, rootward first, so each consecutive pair is one
# bridge over an equal share of the edge.
getMapFromStateSampleBridge <- function(state_sample, edge_length, uni){
  maps <- vector("list", length(state_sample))
  for(edge_i in seq_along(state_sample)){
    path <- state_sample[[edge_i]]
    n_sub <- length(path) - 1L
    sub_length <- edge_length[edge_i] / n_sub
    segments <- vector("list", n_sub)
    for(j in seq_len(n_sub)){
      drawn <- sampleBranchBridge(path[j], path[j + 1L], sub_length, uni)
      if(is.null(drawn)) return(NULL)
      segments[[j]] <- drawn
    }
    edge_map <- unlist(segments, use.names = TRUE)
    # a shift can fall exactly on the join between two sub-edges, so the runs are
    # collapsed once more across the whole edge
    runs <- rle(as.numeric(names(edge_map)))
    grouping <- rep(seq_along(runs$lengths), runs$lengths)
    maps[[edge_i]] <- stats::setNames(
      as.numeric(tapply(edge_map, grouping, sum)), runs$values)
  }
  maps
}
