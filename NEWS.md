# OUwie 3.0.4

* `hOUwie()` now starts sigma^2 from the mean squared independent contrast, as
  `OUwie()` does, instead of the raw trait variance. On tall trees the old start
  could drive alpha to its upper bound and stall far below BM1.
* `hOUwie()` default `root.p` is now `"maddfitz"` (was `"yang"`).
* Fixed the OU expected means and variances along paths with more than one
  alpha, which feed the model-averaged values from `getModelAvgParams()`.
* `hOUwie.walk()` no longer fails on models without alpha (BM1, BMV).

# OUwie 3.0.3

* `hOUwie()` likelihoods changed for all models (importance sampling replaces
  the old subset sum; cache bug fixed). Results from 2.x are not comparable
  and should be refit.

# OUwie 3.01

* Fixed the multiple-alpha (OUMA, OUMVA) calculation. See the Calculation
  Update vignette.
