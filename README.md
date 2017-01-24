# lcopula
## Liouville copulas

Various functions to estimate one of five Liouville copula families as defined in the `copula` package. 
Since the models only admit closed-form survival functions for the integer-valued cased, the optimization routine proceeds
pointwise over grid values for the parameter `alpha` and optimizes over `rho`. The procedure is computationally intensive.


### Possible functions of interest include

- `rliouv` to generate from the Liouville copula, but also `dliouv`, `pliouv` for density and distribution function
- `dliouvm`, `sliouvm`, `isliouvm`, `pliouvm` for marginals of the Liouville copulas
- `liouv.maxim` to fit a copula model to data for some possible combinations of `alpha`. It is possible to rely on a Monte-Carlo approximation of the inverse survival function to speed up computations. This is the default and is heavily recommended
- `liouv.maxim.mm` for method-of-moments based estimation: for the Gumbel and Clayton families, explicit formulas for Kendall's tau can be used to define a moment-based estimator
- `pickands.liouv` and `pickands.plot` for the Pickands dependence function of extremal attractors of Liouville copulas (respectively the negative scaled Dirichlet when `CDA="C"` and positive scaled Dirichlet when `CDA="S"`)
- `hbvevdliouv` for the spectral densities of extremal attractors
- `liouv.Tau` for the theoretical Kendall's tau value for a given set of parameters from the Gumbel or Clayton Liouville copulas
- `K.plot` for Kendall dependence plots

### Datasets

Three datasets, `airquality`, `danube` and `nutrient` are included with the package.

### Installation 
The package is available from CRAN. For up-to-date versions, you can install it from from Github via

```R
devtools::install_github("lbelzile/lcopula")
```

after installing `devtools`.
