# hOUwie optimizer comparison

## Question

Does changing the derivative-free optimizer improve hOUwie's joint discrete/continuous optimization once common random numbers (CRN) make comparisons between nearby parameter vectors less noisy?

## Design

The benchmark simulates a 96-tip tree with four latent states (`1A`, `2A`, `1B`, and `2B`) but exposes only observed states `1` and `2` to the fit. The realized history has 23 changes, observed tip counts of 64/32, and latent-category tip counts of 67/29. Observed-state transition rates differ between latent categories (0.8 versus 1.5), latent transitions are asymmetric (A-to-B 0.3 and B-to-A 0.2), and the two latent categories have subtly different OU optima (2.5 versus 3.7).

Every optimizer used:

- the same four starting vectors and the same CRN seed within each paired start;
- `nSim = 25`, `sample_tips = FALSE`, `sample_nodes = FALSE`, and `adaptive_sampling = FALSE`;
- parameter-dependent stochastic maps at every objective evaluation (CRN fixes the underlying uniform draws, not the maps);
- native convergence criteria rather than an evaluation-count stopping rule.

SBPLX and COBYLA used NLopt. NOMAD used MADS through the CRAN `crs` package. NOMAD's required evaluation ceiling was set to 1,000,000 as a nonbinding safeguard; all runs stopped after 513--823 evaluations because the mesh converged. SBPLX and COBYLA stopped on relative function tolerance.

Each completed fit was evaluated on 30 common, fresh random-number streams at `nSim = 25`. It was then evaluated without further optimization on 10 fresh streams at `nSim = 250`.

## Results

### Fresh-stream validation at the optimization setting (`nSim = 25`)

| Optimizer | Mean log likelihood | Between-start SD | Mean evaluations | Mean time (s) |
|---|---:|---:|---:|---:|
| NOMAD | -52.149 | 0.435 | 710.0 | 212.3 |
| SBPLX | -52.551 | 0.695 | 545.5 | 162.2 |
| COBYLA | -57.786 | 4.135 | 158.0 | 47.2 |

NOMAD improved the four-start mean by 0.402 log-likelihood units over SBPLX and reduced between-start variation, at a cost of about 30% more evaluations and wall time. NOMAD beat its paired SBPLX run for three of four starts; the paired differences were -0.086, +0.786, +0.905, and +0.003.

### Higher-coverage validation (`nSim = 250`)

| Optimizer | Mean log likelihood | Between-start SD |
|---|---:|---:|
| NOMAD | -49.980 | 0.266 |
| SBPLX | -50.338 | 0.522 |
| COBYLA | -55.068 | 3.714 |

The ordering survives the tenfold increase in map count. NOMAD's paired advantages over SBPLX were -0.085, +0.709, +0.806, and +0.003, for a mean gain of 0.358. However, the best NOMAD solution (-49.705) and best SBPLX solution (-49.716) differ by only 0.012. NOMAD's benefit is therefore mostly better reliability across starts, not a higher best-achievable likelihood.

## Parameter warning

Six of the twelve fits had at least one parameter within 2% of a bound. Both NOMAD and SBPLX found solutions with discrete rates at their lower or upper bounds, and some fits effectively eliminated a latent category. The correctly ordered generating vector also scored substantially below the fitted vectors (-63.473 at `nSim = 250`). A single simulated realization need not have its MLE at the generating values, but the size of this gap together with the boundary solutions indicates a ridge or weak-identifiability problem. Increasing `nSim` reduced validation noise but did not remove it.

## Recommendation

NOMAD is a useful optional optimizer, especially for a single start: it was more consistent and modestly better on average. It is not yet compelling as the default because it is about 30% slower, adds a compiled external dependency, and offers essentially no improvement over the best of four SBPLX starts. COBYLA should not be added based on this test.

Before changing the public API, the stronger test is a replicate-level recovery benchmark across multiple simulated trees and datasets. That would separate optimizer robustness from this particular realization's identifiability and reveal whether NOMAD improves parameter recovery, not only likelihood.

## Reproduction

- Main benchmark: `tools/benchmarks/optimizer-comparison.R`
- High-coverage validation: `tools/benchmarks/optimizer-high-nsim-validation.R`
- Moderate-case results: `tools/benchmarks/optimizer-comparison-results-moderate/`
- Earlier extreme stress-test results: `tools/benchmarks/optimizer-comparison-results/`

The earlier stress test produced 96 realized history changes and much stronger finite-map pathologies; it is retained as an extreme stress case rather than used for the recommendation above.
