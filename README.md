
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CURobustML

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/youmisuk/CURobustML.svg?branch=master)](https://travis-ci.com/youmisuk/CURobustML)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/youmisuk/CURobustML?branch=master&svg=true)](https://ci.appveyor.com/project/youmisuk/CURobustML)
[![Codecov test
coverage](https://codecov.io/gh/youmisuk/CURobustML/branch/master/graph/badge.svg)](https://codecov.io/gh/youmisuk/CURobustML?branch=master)
<!-- badges: end -->

A family of machine learning estimators for handling unmeasured
cluster-level confounders in R. It provides a general approach to
estimate causal effects in the presence of cluster-level unmeasured
confounders in multilevel observational data. In particular, we leverage
modern ML methods and exploit a fundamental nature regarding
cluster-level unmeasured confounders to estimate the conditional average
treatment effect (CATE) and the average treatment effect (ATE).

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
```

Or, use the setup function.

``` r
CUClusterML.setup()
```

## Example

  - The doubly robust estimator with proxy regression with two-level
    data

<!-- end list -->

``` r
library(CURobustML)

DRcomb.rslt <- DRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                      X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)

summary(DRcomb.rslt)

# with final predictions
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRcomb.rslt$Z.hat,
   Y1.hat=DRcomb.rslt$Y1.hat, Y0.hat=DRcomb.rslt$Y0.hat, data=twolevel_data)

# with predictions from ensemble glm
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRcomb.rslt$Z.hats$Ztest.hat_glm,
   Y1.hat=DRcomb.rslt$Y1.hats$Y1test.hat_glm, Y0.hat=DRcomb.rslt$Y0.hats$Y0test.hat_glm, data=twolevel_data)

# with predictions from deep learning
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRcomb.rslt$Z.hats$Ztest.hat_dl,
   Y1.hat=DRcomb.rslt$Y1.hats$Y1test.hat_dl, Y0.hat=DRcomb.rslt$Y0.hats$Y0test.hat_dl, data=twolevel_data)
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
