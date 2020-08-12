
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CURobustML

<!-- badges: start -->

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

## Example

  - The doubly robust estimator with proxy regression with two-level
    data

<!-- end list -->

``` r
library(CURobustML)
#> Registered S3 methods overwritten by 'lme4':
#>   method                          from
#>   cooks.distance.influence.merMod car 
#>   influence.merMod                car 
#>   dfbeta.influence.merMod         car 
#>   dfbetas.influence.merMod        car

DRcomb.rslt <- DRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                      X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)
#> [1] "No crossFitting declared. Overriding mCrossfit and K variables to 1 both."
#> [1] "1 round out of  1  and estimated values: 2.001,0.974"
#> [1] "1 round out of  1  and estimated SE: 0.071,0.07"

summary(DRcomb.rslt)
#> Coefficients: 
#>  
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.00097    0.07101   28.18   <2e-16 ***
#> W1           0.97421    0.07008   13.90   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# with final predictions
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRcomb.rslt$Z.hat,
   Y1.hat=DRcomb.rslt$Y1.hat, Y0.hat=DRcomb.rslt$Y0.hat, data=twolevel_data)
#>              Estimate Std. Error
#> (Intercept) 2.0009661 0.07100601
#> W1          0.9742097 0.07008262

# with predictions from ensemble glm
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRcomb.rslt$Z.hats$Ztest.hat_glm,
   Y1.hat=DRcomb.rslt$Y1.hats$Y1test.hat_glm, Y0.hat=DRcomb.rslt$Y0.hats$Y0test.hat_glm, data=twolevel_data)
#>              Estimate Std. Error
#> (Intercept) 2.0007282 0.07068592
#> W1          0.9707472 0.06976670

# with predictions from deep learning
DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRcomb.rslt$Z.hats$Ztest.hat_dl,
   Y1.hat=DRcomb.rslt$Y1.hats$Y1test.hat_dl, Y0.hat=DRcomb.rslt$Y0.hats$Y0test.hat_dl, data=twolevel_data)
#>             Estimate Std. Error
#> (Intercept) 1.940238  0.2720784
#> W1          1.637338  0.2685402
```

  - The double demeaning estimator with cross-classified data

<!-- end list -->

``` r
DDcomb.rslt <- DDcomb(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1),
                      X=crossclassified_data[, c("X1", "X2", "X3", "W1", "Q1")],
                      ID=crossclassified_data$f12id, data=crossclassified_data)
#> [1] "No crossFitting declared. Overriding mCrossfit and K variables to 1 both."
#> [1] "1 round out of  1  and estimated values: 1.961,0.979"
#> [1] "1 round out of  1  and estimated SE: 0.055,0.04"

summary(DDcomb.rslt)
#> Coefficients: 
#>  
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.96075    0.05533   35.44   <2e-16 ***
#> W1           0.97940    0.04027   24.32   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# with final predictions
DD(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1), ID=crossclassified_data$f12id,
   Z.hat=DDcomb.rslt$Z.hat, Y0.hat=DDcomb.rslt$Y0.hat, data=crossclassified_data)
#>              Estimate Std. Error
#> (Intercept) 1.9607486  0.0553275
#> W1          0.9794042  0.0402655

# with predictions from glm
DD(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1), ID=crossclassified_data$f12id,
   Z.hat=DDcomb.rslt$Z.hats$Ztest.hat_glm, Y0.hat=DDcomb.rslt$Y0.hats$Y0test.hat_glm, data=crossclassified_data)
#>              Estimate Std. Error
#> (Intercept) 1.9746300 0.04938442
#> W1          0.9976662 0.03594905

# with predictions from deep learning
DD(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1), ID=crossclassified_data$f12id,
   Z.hat=DDcomb.rslt$Z.hats$Ztest.hat_dl, Y0.hat=DDcomb.rslt$Y0.hats$Y0test.hat_dl, data=crossclassified_data)
#>              Estimate Std. Error
#> (Intercept) 1.9607486  0.0553275
#> W1          0.9794042  0.0402655
```
