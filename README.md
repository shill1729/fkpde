
# fkpde

<!-- badges: start -->
<!-- badges: end -->

This package provides finite difference solvers for PDEs in Feynman-Kac problems concerning stochastic diffusions.


## Installation

You can install the current GitHub version via devtools

``` r
devtools::install_github("shill1729/fkpde")
```

## Example
Basic Black-Scholes PDE
```r
# Define dynamics
dynamics <- list(function(t, x) 0.05*x,
                 function(t, x) 0.55*x
                 )
# Define conditional expectation
problem <- list(function(t, x) 0,
                function(t, x) 0, 
                function(x) pmax(x-100, 0)
                )
# Define region and control parameters
region <- c(1, 0, 200)
control <- list(N = 100, M = 100, variational = FALSE, output = "grid")
fk <- solvePDE(dynamics, problem, region, control)
plot(fk$solution[control$N+1, ], type = "l")
```

