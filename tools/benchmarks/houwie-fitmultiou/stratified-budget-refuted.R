#!/usr/bin/env Rscript
#
# REFUTED - kept so the idea is not tried again. Result in
# results/cap-surface/stratified-budget-refuted.txt.
#
# Does reserving budget per descendant regime fix the four-regime near-optimum
# error? The package reduction concatenates the pushed messages from all K
# descendant regimes and merges the pile globally, so a mass-weighted merge can
# delete a whole regime's subpopulation. GPB2/IMM never merge across discrete
# states. This gives each descendant regime a floor of budget/K.
#
# It makes the error worse in all sixteen comparisons - roughly two to five times
# worse near the optimum at every cap. The premise was wrong: components from
# different regimes sit at different optima, so Runnalls' cost between them is
# already high and the global merge was not crossing regimes to begin with.
# Forcing a per-regime budget only compels within-regime merges that the global
# merge had correctly declined, so it spends accuracy to buy a guarantee that was
# never being violated.
#
# Running it with stratify = FALSE reproduces houwiePruningLik exactly, which is
# what makes the paired comparison valid.
suppressPackageStartupMessages({library(ape); library(phytools)
  pkgload::load_all("/Users/jboyko/OUwie", quiet=TRUE, export_all=FALSE, helpers=FALSE)})

pruneStratified <- function(phy, tip_state_sets, tip_values, Q, alpha, sigma.sq,
                            theta, root.p, resolution = 1L,
                            max_components = Inf, stratify = TRUE){
  n_tip <- length(phy$tip.label); nStates <- nrow(Q)
  red <- function(m, cap) OUwie:::reduceMessageMixture(m, cap, 0, NULL)
  tcache <- new.env(parent = emptyenv())
  transitionFor <- function(d){ key <- sprintf("%.14g", d)
    if(is.null(tcache[[key]])) tcache[[key]] <- as.matrix(expm::expm(Q*d))
    tcache[[key]] }
  sendUp <- function(messages, edge_length, allowed){
    interval <- edge_length/resolution; transition <- transitionFor(interval)
    half <- interval/2
    composed <- lapply(seq_len(nStates), function(a) lapply(seq_len(nStates),
      function(d) OUwie:::composeIntervalPair(a, d, half, alpha, sigma.sq, theta)))
    for(step in seq_len(resolution)){
      updated <- vector("list", nStates)
      for(a in seq_len(nStates)){
        pieces <- list()
        for(d in seq_len(nStates)){
          if(step == 1L && !is.null(allowed) && !(d %in% allowed)) next
          ls <- log(transition[a, d]); if(!is.finite(ls)) next
          mo <- composed[[a]][[d]]
          p <- OUwie:::pushMessage(messages[[d]], mo[["A"]], mo[["B"]], mo[["V"]])
          p$g <- p$g + ls
          pieces[[length(pieces)+1L]] <- p
        }
        updated[[a]] <- if(!length(pieces)) OUwie:::emptyMessage() else {
          if(stratify && is.finite(max_components) && length(pieces) > 1L){
            share <- max(1L, as.integer(ceiling(max_components/length(pieces))))
            pieces <- lapply(pieces, red, cap = share)
          }
          red(OUwie:::bindMessages(pieces), max_components)
        }
      }
      messages <- updated
    }
    messages
  }
  node <- vector("list", n_tip + phy$Nnode)
  for(anc in unique(phy$edge[,1])){
    combined <- NULL
    for(e in which(phy$edge[,1] == anc)){
      des <- phy$edge[e,2]
      if(des <= n_tip){
        obs <- list(g=0, h=unname(tip_values[des]), k=Inf)
        child <- rep(list(obs), nStates); allowed <- tip_state_sets[[des]]
      } else { child <- node[[des]]; allowed <- NULL }
      up <- sendUp(child, phy$edge.length[e], allowed)
      combined <- if(is.null(combined)) up else lapply(seq_len(nStates),
        function(s) red(OUwie:::multiplyMessages(combined[[s]], up[[s]]), max_components))
    }
    node[[anc]] <- combined
  }
  rm_ <- node[[n_tip+1L]]
  per <- vapply(seq_len(nStates), function(s){ m <- rm_[[s]]; st <- theta[s]
    OUwie:::logSumExpSafe(m$g + st*(m$h - 0.5*m$k*st)) }, numeric(1))
  OUwie:::logSumExpSafe(log(root.p/sum(root.p)) + per)
}

set.seed(5208); phy <- pbtree(n=8, scale=10)
Q2 <- matrix(c(-.3,.3,.3,-.3),2,2,byrow=TRUE)
sim <- hOUwie.sim(phy, Q=Q2, root.freqs=c(1,0), alpha=c(.5,.5), sigma.sq=c(.2,.2), theta0=5, theta=c(5,10))
phy <- reorder.phylo(phy,"pruningwise")
values <- setNames(sim$data[,3], sim$data[,1])[phy$tip.label]
obs <- setNames(as.integer(sim$data[,2]), sim$data[,1])[phy$tip.label]
amb <- lapply(obs, function(s) if(s==1L) c(1L,2L) else c(3L,4L)); names(amb) <- names(obs)
Q4 <- matrix(0,4,4); Q4[1,2]<-Q4[2,1]<-Q4[3,4]<-Q4[4,3]<-.15
Q4[1,3]<-Q4[3,1]<-Q4[2,4]<-Q4[4,2]<-.30; diag(Q4) <- -rowSums(Q4)

scans <- list(
 "hidden alpha ratio"=list(axis=exp(seq(log(1),log(20),length.out=13)),
   build=function(v) list(a=c(.25,.25*v,.25,.25*v), s=rep(.3,4), t=c(5,5,12,12))),
 "hidden theta separation"=list(axis=seq(0,14,length.out=13),
   build=function(v) list(a=c(.25,1,.25,1), s=c(.15,.6,.15,.6), t=c(5,5,5+v,5+v))))

for(nm in names(scans)){
  sc <- scans[[nm]]
  ev <- function(cap, strat) vapply(sc$axis, function(v){ p <- sc$build(v)
    pruneStratified(phy, amb, values, Q4, p$a, p$s, p$t, rep(.25,4), 1L, cap, strat)}, numeric(1))
  exact <- ev(Inf, FALSE); below <- max(exact)-exact
  cat(sprintf("\n=== %s (8 tips, K=4) ===\n", nm))
  for(cap in c(8,16,32,64)){
    for(strat in c(FALSE, TRUE)){
      e <- ev(cap, strat) - exact
      cat(sprintf("  cap=%-3d stratify=%-5s  within2=%9.3e  within10=%9.3e  all=%8.3f\n",
        cap, strat, max(abs(e[below<=2])), max(abs(e[below<=10])), max(abs(e))))
    }
  }
}
