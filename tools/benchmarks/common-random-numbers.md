# Common-random-number benchmark (2026-07-31)

Question: does resetting the underlying random-number stream for each likelihood
evaluation improve hOUwie optimization without fixing stochastic maps?

In all comparisons, maps were regenerated for the current parameters. Only the
random-number stream was shared within an optimizer start. Current and CRN fits
received equal settings and were evaluated afterward on fresh random streams.

## Documented ER + OUM example

Configuration: `tworegime`, `rate.cat = 1`, `nSim = 25`, with
`sample_nodes = FALSE` and `adaptive_sampling = FALSE`.

CRN reduced the variance of 1% local parameter contrasts by approximately 35x
to 22,000x. All eight current and CRN searches converged to the same parameter
vector and the same fresh-stream log-likelihood (`-30.9956`). This demonstrated
smoother comparisons but no change in the fitted result on the easy problem.

## Difficult documented CID model

Configuration: `tworegime`, ER + OUM, `rate.cat = 2`, `null.model = TRUE`,
`nSim = 25`, and default sampling flags. Six paired SBPLX searches were allowed
to converge without an evaluation cutoff, then each fit was evaluated on 30
fresh random streams.

| metric | current | CRN |
|---|---:|---:|
| fresh-stream mean log-likelihood | -27.9244 | -27.2323 |
| fresh-stream median log-likelihood | -28.0759 | -27.1515 |
| between-run SD | 0.6145 | 0.2669 |
| mean evaluations | 524.0 | 438.5 |
| mean elapsed seconds | 75.24 | 61.48 |
| mean optimization optimism | 0.3798 | 0.1286 |

CRN won four paired starts and lost two. Mean and median fresh-stream gains were
0.692 and 0.939 log-likelihood units. It used 16.3% fewer evaluations and 18.3%
less wall time. All current searches stopped on parameter-step tolerance, while
all CRN searches stopped on function tolerance.

Several hidden-state transition rates still approached their bounds and
parameter estimates varied across starts. CRN improved this search but did not
remove identifiability or local-optimum concerns. This supports an experimental
opt-in argument, not a default change.
