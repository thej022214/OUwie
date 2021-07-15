[![Build Status](https://travis-ci.org/thej022214/OUwie.svg)](https://travis-ci.org/thej022214/OUwie)

OUwie is an R package for using Brownian motion and Ornstein-Uhlenbeck models for trait evolution.
Its friendly webpage is at http://thej022214.github.io/OUwie/; its source code is at https://github.com/thej022214/OUwie/. 

Some of the features:

* Brownian motion models that allow the rate (sigma-squared) to vary over the tree
* Ornstein-Uhlenbeck models that allow the rate, optima (theta), and/or strength of pull (alpha) to vary over the tree
* Uncertainty estimation using contour plots to find potential ridges
* Simulation functions
* Automatic testing of some identifiability issues using methods from Ho and Ané (2014)
* Ancestral state estimation under all these models (though use substantial caution)
* Use of measurement error at the tips

Some of its caveats:

* It is univariate (a single trait) only
* For multiple rate models, it requires some mapping of regimes (stochastic character mapping of a discrete state, using node labels for regimes on trees, etc.).
* It warns you about models that are very complex for what your data may allow, but it will let you run them
* Optimization can be a difficult problem -- it tries its best, and will announce failures when it notices them, but still be careful

This is the bleeding edge version: you can install it with `remotes::install_github("thej022214/OUwie")` [install the remotes package from CRAN first]

