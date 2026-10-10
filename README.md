Average Shared Variance-Based Relative Importance Analysis for Fixed-Effects Regression Models

`feols.asv` is an R function for evaluating the relative importance of correlated explanatory variables in fixed-effects regression models fitted using the `feols()` function from the `fixest` package. Based on the **Average Shared Variance (ASV)** framework, it decomposes the within-model \(R^2\) into unique and shared components and quantifies the individual contribution of each predictor by equally allocating shared explained variation among the predictors involved.

The function supports models with multiple dimensions of fixed effects, such as county and year fixed effects, while retaining the specified fixed-effects structure across predictor subsets. It provides an order-independent and additive assessment of predictor importance, helping researchers evaluate the relative contributions of correlated explanatory variables in fixed-effects regression analyses.

## Installation

Install the package from GitHub using:

```r
remotes::install_github("laijiangshan/feols.asv")
```

## Related R packages

The ASV framework is relevant to a broader ecosystem of R packages developed for relative importance analysis and explained-variation partitioning in ecological and environmental research, including:

- `rdacca.hp` — relative importance analysis and variation partitioning in canonical analyses
- `glmm.hp` — relative importance analysis for generalized linear mixed models
- `gam.hp` — relative importance analysis for generalized additive models
- `phylolm.hp` — relative importance analysis for phylogenetic linear models
- `spatialreg.hp` — relative importance analysis for spatial autoregressive models
- `asv.hp` — ASV-based relative importance analysis and commonality analysis

## License and Contact

This package is distributed under the license specified in the package repository.

For questions or further information, please contact [lai@njfu.edu.cn](mailto:lai@njfu.edu.cn).
