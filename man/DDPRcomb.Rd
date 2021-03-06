% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combFun.R
\name{DDPRcomb}
\alias{DDPRcomb}
\title{Double demeaning estimator with proxy regressions with treatment and outcome learners}
\usage{
DDPRcomb(
  Y,
  Z,
  X,
  interZ = formula(~1),
  ID,
  data,
  library = c("glm", "deeplearning"),
  crossFitting = FALSE,
  K = 5,
  mCrossFit = 100
)
}
\arguments{
\item{Y}{continuous outcome variable}

\item{Z}{binary treatment indicator, 1 - treatment, 0 - control}

\item{X}{vector, matrix, or dataframe containing measured confounders}

\item{interZ}{formula that contains the variables that "interact"  with the treatment. "1" will be always added. The default is no interaction, i.e., formula = formula(~1).}

\item{ID}{cluster identifier}

\item{data}{dataframe containing the variables in the model}

\item{library}{character vector of prediction algorithms. The available methods are  glm, deeplearning, gbm, and randomForests. The default methods are glm and deeplearning. Three types of glm include: glm with cluster dummies, glm with cluster-constant components of individual-level covariates, and glmm with random effects of clusters.}

\item{crossFitting}{whether to do cross fitting. The default is FALSE, and currently it is not available.}

\item{K}{number of folds. The default is 5, and currently it is not available.}

\item{mCrossFit}{number of cross fitting. The default is FALSE, and currently it is not available.}
}
\value{
An \code{DDPRcomb} with the following elements:
   \item{coef.ER}{vector of the coefficients for prediction algorithms in the treatment model}
   \item{coef.OR}{vector of the coefficients for prediction algorithms in the outcome model}
   \item{Estimate}{estimates and standard errors of treatment effects}
   \item{Z.hat}{final weighted prediction for the treatment}
   \item{Y1.hat}{final weighted prediction for the outcome among treated units}
   \item{Y0.hat}{final weighted prediction for the outcome among control units}
   \item{Z.hats}{all the predictions for the treatment from prediction algorithms}
   \item{Y1.hats}{all the predictions for the outcome among treated units from prediction algorithms}
   \item{Y0.hats}{all the predictions for the outcome among control units from prediction algorithms}
}
\description{
Double demeaning estimator with proxy regressions with treatment and outcome learners
}
\examples{
# two-level data
DDPRcomb.rlst <- DDPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
 X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id,
 library="glm", data=twolevel_data)
summary(DDPRcomb.rlst)

\donttest{
# two-level data with default algorithms
DDPRcomb.rlst2 <- DDPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
 X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id,
 library=c("glm", "deeplearning"), data=twolevel_data)
summary(DDPRcomb.rlst2)

# cross-classified data with default algorithms
DDPRcomb.rlst3 <- DDPRcomb(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1),
 X=crossclassified_data[, c("X1", "X2", "X3", "W1", "Q1")], ID=crossclassified_data$f12id,
 library=c("glm", "deeplearning"), data=crossclassified_data)
summary(DDPRcomb.rlst3)
}
}
