
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CURobustML

A family of machine learning estimators for handling unmeasured
cluster-level confounders in R. It provides a general approach to
estimate causal effects in the presence of cluster-level unmeasured
confounders in multilevel observational data. In particular, we leverage
modern ML methods and exploit a fundamental nature regarding
cluster-level unmeasured confounders to estimate the conditional average
treatment effect (CATE) and the average treatment effect (ATE). See Suk
and Kang (2020) \<doi:
[10.31234/osf.io/t7vbz](https://doi.org/10.31234/osf.io/t7vbz)\> for
details.

## Installation

You can install the released version of CURobustML from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CURobustML")
```

You should also start and connect to H2O instance via package h2o to
implement CURobustML.

``` r
library(h2o)
localH2O = h2o.init()
h2o.removeAll() # clean slate - just in case the cluster was already running
```

Or, use the setup function.

``` r
CURobustML.setup()
```

## Example

  - The doubly robust estimator with proxy regression with two-level
    data

<!-- end list -->

``` r
library(CURobustML)

DRPRcomb.rslt <- DRPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                      X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)

summary(DRPRcomb.rslt)

# with final predictions
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRPRcomb.rslt$Z.hat,
   Y1.hat=DRPRcomb.rslt$Y1.hat, Y0.hat=DRPRcomb.rslt$Y0.hat, data=twolevel_data)

# with predictions from glm wtih fixed effects of clusters
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRPRcomb.rslt$Z.hats$Ztest.hat_glm_1,
   Y1.hat=DRPRcomb.rslt$Y1.hats$Y1test.hat_glm_1, Y0.hat=DRPRcomb.rslt$Y0.hats$Y0test.hat_glm_1, data=twolevel_data)

# with predictions from glm wtih random effects of clusters
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRPRcomb.rslt$Z.hats$Ztest.hat_glm_3,
   Y1.hat=DRPRcomb.rslt$Y1.hats$Y1test.hat_glm_3, Y0.hat=DRPRcomb.rslt$Y0.hats$Y0test.hat_glm_3, data=twolevel_data)
```

  - The double demeaning estimator with cross-classified data

<!-- end list -->

``` r
DDcomb.rslt <- DDcomb(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1),
                      X=crossclassified_data[, c("X1", "X2", "X3", "W1", "Q1")],
                      ID=crossclassified_data$f12id, data=crossclassified_data)

summary(DDcomb.rslt)

# with final predictions
DD(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1), ID=crossclassified_data$f12id,
   Z.hat=DDcomb.rslt$Z.hat, Y0.hat=DDcomb.rslt$Y0.hat, data=crossclassified_data)

# with predictions from glm
DD(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1), ID=crossclassified_data$f12id,
   Z.hat=DDcomb.rslt$Z.hats$Ztest.hat_glm, Y0.hat=DDcomb.rslt$Y0.hats$Y0test.hat_glm, data=crossclassified_data)

# with predictions from deep learning
DD(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1), ID=crossclassified_data$f12id,
   Z.hat=DDcomb.rslt$Z.hats$Ztest.hat_dl, Y0.hat=DDcomb.rslt$Y0.hats$Y0test.hat_dl, data=crossclassified_data)
```
