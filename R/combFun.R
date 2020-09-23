
#' Doubly robust estimator with proxy regressions with treatment and outcome learners
#'
#' @param Y continuous outcome variable
#' @param Z binary treatment indicator, 1 - treatment, 0 - control
#' @param X vector, matrix, or dataframe containing measured confounders
#' @param interZ formula that contains the variables that "interact"  with the treatment. "1" will be always added. The default is no interaction, i.e., formula = formula(~1).
#' @param ID cluster identifier
#' @param data dataframe containing the variables in the model
#' @param library character vector of prediction algorithms. The available methods are  glm, deeplearning, gbm, and randomForests. The default methods are glm and deeplearning. Three types of glm include: glm with cluster dummies, glm with cluster-constant components of individual-level covariates, and glmm with random effects of clusters.
#' @param crossFitting whether to do cross fitting. The default is FALSE, and currently it is not available.
#' @param K number of folds. The default is 5, and currently it is not available.
#' @param mCrossFit number of cross fitting. The default is FALSE, and currently it is not available.
#'
#' @return
#' An \code{DRPRcomb} with the following elements:
#'    \item{coef.ER}{vector of the coefficients for prediction algorithms in the treatment model}
#'    \item{coef.OR}{vector of the coefficients for prediction algorithms in the outcome model}
#'    \item{Estimate}{estimates and standard errors of treatment effects}
#'    \item{Z.hat}{final weighted prediction for the treatment}
#'    \item{Y1.hat}{final weighted prediction for the outcome among treated units}
#'    \item{Y0.hat}{final weighted prediction for the outcome among control units}
#'    \item{Z.hats}{all the predictions for the treatment from prediction algorithms}
#'    \item{Y1.hats}{all the predictions for the outcome among treated units from prediction algorithms}
#'    \item{Y0.hats}{all the predictions for the outcome among control units from prediction algorithms}
#' @export
#'
#' @importFrom stats formula median model.matrix binomial predict
#' @importFrom nnls nnls
#' @importFrom caret createFolds
#' @importFrom WeightIt make_full_rank
#' @importFrom lme4 glmer lmer glmerControl lmerControl
#' @importFrom h2o as.h2o h2o.glm h2o.gbm h2o.deeplearning h2o.randomForest
#'
#' @examples
#' # with two-level data
#' DRPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
#'  X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)
#'
#' \donttest{
#' # with cross-classified data
#' DRPRcomb(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1),
#'  X=crossclassified_data[, c("X1", "X2", "X3", "W1", "Q1")], ID=crossclassified_data$f12id,
#'  data=crossclassified_data)
#' }
DRPRcomb <- function(Y, Z, X, interZ=formula(~1), ID, data, library=c("glm", "deeplearning"), crossFitting=FALSE,
                     K = 5, mCrossFit = 100) {

  if(length(library) != sum(library %in% c("glm", "deeplearning", "gbm", "randomForests"))) {
    stop("library is not in the appropriate format. Choose methods among glm, gbm, randomForests, and deeplearning.")
  }

  if (crossFitting) {
    stop("Cross fitting has not yet been implemented.")
  }

  if(!crossFitting) {
    print("No crossFitting declared. Overriding mCrossfit and K variables to 1 both.")
    mCrossFit = 1; K = 1
  }

  # variables
  IDfactor = as.factor(ID)
  IDnum = as.numeric(IDfactor)
  Zfactor = as.factor(Z)
  Z = as.numeric(as.character(Z))

  lvl2var <- colSums(apply(as.matrix(X), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1 <- as.matrix(X)[, !lvl2var]
  W_lvl2 <- as.matrix(X)[, lvl2var]
  X_lvl1_num <- as.matrix(as.matrix(X_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
  X_lvl1_mean <- apply(as.matrix(X)[, !lvl2var], 2, function(y) ave(y, IDnum))

  X_lvl1_num_sq <- X_lvl1_num^2
  colnames(X_lvl1_num_sq) <- paste0(colnames(X_lvl1_num_sq), "_sq")

  # datasets for the propensity scores
  h2oDataE= as.h2o(data.frame(Zfactor, X_lvl1, W_lvl2, IDfactor, X_lvl1_num_sq))
  h2oDataE2 = as.h2o(data.frame(Zfactor, X_lvl1, W_lvl2, X_lvl1_mean, X_lvl1_num_sq))

  VarAll = data.frame(X_lvl1,W_lvl2)
  VarInt.m <- model.matrix(~ .^2, data=VarAll)
  VarInt.m2 <- WeightIt::make_full_rank(VarInt.m, with.intercept = TRUE)
  REdataE = data.frame(VarInt.m2, X_lvl1_num_sq, Zfactor, IDfactor)

  covFrml <- paste(colnames(VarInt.m2), collapse="+")
  glmerFrml <- formula(paste("Zfactor ~", covFrml, " + (1 | IDfactor)"))

  # datasets for the outcome
  h2oAll = as.h2o(data.frame(Y, Z, X_lvl1, W_lvl2 ,IDfactor, X_lvl1_num_sq))
  h2oAll_Z1 = h2oAll_Z0 = h2oAll; h2oAll_Z1$Z = 1; h2oAll_Z0$Z = 0

  h2oAll2 = as.h2o(data.frame(Y, Z, X_lvl1, W_lvl2, X_lvl1_mean, X_lvl1_num_sq))
  h2oAll2_Z1 = h2oAll2_Z0 = h2oAll2; h2oAll2_Z1$Z = 1; h2oAll2_Z0$Z = 0

  VarTmp = data.frame(Z,X_lvl1,W_lvl2)
  VarAll2 = rbind(VarTmp,VarTmp,VarTmp)
  VarAll2$Z[(nrow(VarTmp) + 1):(2*nrow(VarTmp))] = 1
  VarAll2$Z[(2*nrow(VarTmp) + 1):(3*nrow(VarTmp))] = 0
  VarInt2.m <- model.matrix(~ .^2, data=VarAll2)
  VarInt2.m2 <- WeightIt::make_full_rank(VarInt2.m, with.intercept = TRUE)

  Var2Tmp = data.frame(Y, X_lvl1_num_sq, IDfactor)
  REdata = data.frame(VarInt2.m2[1:nrow(VarTmp), ], Var2Tmp)
  REdata_Z1 = data.frame(VarInt2.m2[(nrow(VarTmp) + 1):(2*nrow(VarTmp)), ], Var2Tmp)
  REdata_Z0 = data.frame(VarInt2.m2[(2*nrow(VarTmp) + 1):(3*nrow(VarTmp)), ], Var2Tmp)

  covFrml2 <- paste(colnames(VarInt2.m2), collapse="+")
  lmerFrml <- formula(paste("Y ~", covFrml2, " + (1 | IDfactor)"))

  interZ.mat <- model.matrix(interZ, data=data)

  if(is.null(interZ)) {
    overallEst = overallSE = matrix(0,mCrossFit,1)
  } else {
    overallEst = overallSE = matrix(0,mCrossFit,ncol(interZ.mat))
  }


  for(m in 1:mCrossFit) {

    dataALL = data.frame(matrix(0, nrow(data), ncol(data)))
    colnames(dataALL) = colnames(data)
    Yall = rep(0,length(Y))
    Y1hatall = rep(0,length(Y))
    Y0hatall = rep(0,length(Y))
    Zall = rep(0,length(Z))
    Zhatall = rep(0,length(Z))

    splitIndex = createFolds(IDfactor,k = K,list=FALSE)

    if(is.null(interZ)) {
      overallFold = matrix(0,K,1)
    } else {
      overallFold = matrix(0,K,ncol(interZ.mat))
    }


    for(i in 1:K) {
      trainIndex = which(splitIndex == i)
      if(!crossFitting) {
        testIndex = trainIndex
      } else {
        testIndex = which(splitIndex != i)
      }

      # ::: propensity score model
      h2oDataETrain = h2oDataE[trainIndex,]; h2oDataETest = h2oDataE[testIndex,]
      h2oDataE2Train = h2oDataE2[trainIndex,]; h2oDataE2Test = h2oDataE2[testIndex,]

      h2oAllTrain = h2oAll[trainIndex,]; h2oAllTest_Z1 = h2oAll_Z1[testIndex,]; h2oAllTest_Z0 = h2oAll_Z0[testIndex,]
      h2oAll2Train = h2oAll2[trainIndex,]; h2oAll2Test_Z1 = h2oAll2_Z1[testIndex,]; h2oAll2Test_Z0 = h2oAll2_Z0[testIndex,]

      REdataETrain = REdataE[trainIndex,]; REdataETest = REdataE[testIndex,]
      REdataTrain = REdata[trainIndex,]; REdataTest_Z1 = REdata_Z1[testIndex,]; REdataTest_Z0 = REdata_Z0[testIndex,]

      fitER_glm <- fitER_dl <- fitER_gbm <- fitER_rf <- Ztest.hat_glm <- Ztest.hat_dl <- Ztest.hat_gbm <- Ztest.hat_rf <-
        Ztrain.hat_glm <- Ztrain.hat_dl <- Ztrain.hat_gbm <- Ztrain.hat_rf <- NULL

      if ("glm" %in% library) {
        # ensemble glm for propensity score training
        if(length(2:(ncol(h2oDataETrain) -ncol(X_lvl1_num)-1))>=2) {
          fitER_glm_1 = h2o.glm(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,alpha=1,lambda=0,family="binomial",interactions=2:(ncol(h2oDataETrain) -ncol(X_lvl1_num)-1))
          fitER_glm_2 = h2o.glm(x=2:ncol(h2oDataE2Train),y=1,training_frame=h2oDataE2Train,alpha=1,lambda=0,family="binomial",interactions=2:(ncol(h2oDataE2Train) -ncol(X_lvl1_num)))
        } else {
          fitER_glm_1 = h2o.glm(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,alpha=1,lambda=0,family="binomial")
          fitER_glm_2 = h2o.glm(x=2:ncol(h2oDataE2Train),y=1,training_frame=h2oDataE2Train,alpha=1,lambda=0,family="binomial")
        }

        fitER_glm_3 = glmer(glmerFrml, family=binomial(), data=REdataETrain, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5), calc.derivs = FALSE), nAGQ = 0)

        Ztrain.hat_glm_1 = as.numeric(as.data.frame(predict(fitER_glm_1,h2oDataETrain)[,"p1"])$p1)
        Ztrain.hat_glm_2 = as.numeric(as.data.frame(predict(fitER_glm_2,h2oDataE2Train)[,"p1"])$p1)
        Ztrain.hat_glm_3 = as.numeric(predict(fitER_glm_3, type="response"))

        Ztest.hat_glm_1 = as.numeric(as.data.frame(predict(fitER_glm_1,h2oDataETest)[,"p1"])$p1)
        Ztest.hat_glm_2 = as.numeric(as.data.frame(predict(fitER_glm_2,h2oDataE2Test)[,"p1"])$p1)
        Ztest.hat_glm_3 = as.numeric(predict(fitER_glm_3, type="response", REdataETest))

      }

      if ("deeplearning" %in% library) {
        fitER_dl = h2o.deeplearning(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,hidden= rep((ncol(h2oDataETrain) -2) + nlevels(IDfactor),2),epochs=100,adaptive_rate=FALSE,rate=10^-5)
        Ztrain.hat_dl = as.numeric(as.data.frame(predict(fitER_dl,h2oDataETrain)[,"p1"])$p1)
        Ztest.hat_dl = as.numeric(as.data.frame(predict(fitER_dl,h2oDataETest)[,"p1"])$p1)
      }

      if ("gbm" %in% library) {
        fitER_gbm = h2o.gbm(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
        Ztrain.hat_gbm = as.numeric(as.data.frame(predict(fitER_gbm,h2oDataETrain)$pred)$predict)
        Ztest.hat_gbm = as.numeric(as.data.frame(predict(fitER_gbm,h2oDataETest)$pred)$predict)
      }

      if ("randomForests" %in% library) {
        fitER_rf = h2o.randomForest(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
        Ztrain.hat_rf = as.numeric(as.data.frame(predict(fitER_rf,h2oDataETrain)$pred)$predict)
        Ztest.hat_rf = as.numeric(as.data.frame(predict(fitER_rf,h2oDataETrain)$pred)$predict)
      }


      # all predictions
      Ztrain.hats <- cbind(Ztrain.hat_glm_1, Ztrain.hat_glm_2, Ztrain.hat_glm_3, Ztrain.hat_dl, Ztrain.hat_gbm, Ztrain.hat_rf)
      Ztest.hats <- cbind(Ztest.hat_glm_1, Ztest.hat_glm_2, Ztest.hat_glm_3, Ztest.hat_dl, Ztest.hat_gbm, Ztest.hat_rf)

      # nnls
      fit.nnlsER <- nnls(A=Ztrain.hats, b=Z)
      initCoefER <- coef(fit.nnlsER)
      initCoefER[is.na(initCoefER)] <- 0.0

      # normalize so sum(coef) = 1 if possible
      if (sum(initCoefER) > 0) {
        coefER <- initCoefER / sum(initCoefER)
      } else {
        warning("All algorithms have zero weight", call. = FALSE)
        coefER <- initCoefER
      }

      Ztest.comb <- as.numeric(crossprod(t(Ztest.hats[, coefER != 0, drop = FALSE]), coefER[coefER != 0]))

      # ::: outcome model
      fitOR_glm <- fitOR_dl <- fitOR_gbm <- fitOR_rf <- Ytrain.hat_glm <- Ytrain.hat_dl <- Ytrain.hat_gbm <- Ytrain.hat_rf <-
        Y1test.hat_glm <- Y1test.hat_dl <- Y1test.hat_gbm <- Y1test.hat_rf <- Y0test.hat_glm <- Y0test.hat_dl <- Y0test.hat_gbm <- Y0test.hat_rf <- NULL

      if ("glm" %in% library) {
        # ensemble glm for outcome reg training
        if(length(2:(ncol(h2oAllTrain) -ncol(X_lvl1_num)-1))>=2) {
          fitOR_glm_1 = h2o.glm(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAllTrain) -ncol(X_lvl1_num)-1))
          fitOR_glm_2 = h2o.glm(x=2:ncol(h2oAll2Train),y=1,training_frame=h2oAll2Train,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAll2Train) -ncol(X_lvl1_num)))
        } else {
          fitOR_glm_1 = h2o.glm(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain,alpha=1,lambda=0,family="gaussian")
          fitOR_glm_2 = h2o.glm(x=2:ncol(h2oAll2Train),y=1,training_frame=h2oAll2Train,alpha=1,lambda=0,family="gaussian")
        }

        fitOR_glm_3 = lmer(lmerFrml, data=REdataTrain, control=lmerControl(calc.derivs = FALSE))

        Ytrain.hat_glm_1 = as.numeric(as.data.frame(predict(fitOR_glm_1,h2oAllTrain)$pred)$predict)
        Ytrain.hat_glm_2 = as.numeric(as.data.frame(predict(fitOR_glm_2,h2oAll2Train)$pred)$predict)
        Ytrain.hat_glm_3 = as.numeric(predict(fitOR_glm_3))

        Y1test.hat_glm_1 = as.numeric(as.data.frame(predict(fitOR_glm_1,h2oAllTest_Z1)$pred)$predict)
        Y1test.hat_glm_2 = as.numeric(as.data.frame(predict(fitOR_glm_2,h2oAll2Test_Z1)$pred)$predict)
        Y1test.hat_glm_3 = as.numeric(predict(fitOR_glm_3, REdataTest_Z1))

        Y0test.hat_glm_1 = as.numeric(as.data.frame(predict(fitOR_glm_1,h2oAllTest_Z0)$pred)$predict)
        Y0test.hat_glm_2 = as.numeric(as.data.frame(predict(fitOR_glm_2,h2oAll2Test_Z0)$pred)$predict)
        Y0test.hat_glm_3 = as.numeric(predict(fitOR_glm_3, REdataTest_Z0))

      }

      if ("deeplearning" %in% library) {
        fitOR_dl = h2o.deeplearning(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain,hidden= rep((ncol(h2oAllTrain) -2) + nlevels(IDfactor),2), epochs=100,adaptive_rate=FALSE,rate=10^-5)
        Ytrain.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTrain)$pred)$predict)
        Y1test.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTest_Z0)$pred)$predict)
      }

      if ("gbm" %in% library) {
        fitOR_gbm = h2o.gbm(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain, ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
        Ytrain.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTrain)$pred)$predict)
        Y1test.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTest_Z0)$pred)$predict)
      }

      if ("randomForests" %in% library) {
        fitOR_rf = h2o.randomForest(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
        Ytrain.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTrain)$pred)$predict)
        Y1test.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTest_Z0)$pred)$predict)
      }

      # all predictions
      Ytrain.hats <- cbind(Ytrain.hat_glm_1, Ytrain.hat_glm_2, Ytrain.hat_glm_3, Ytrain.hat_dl, Ytrain.hat_gbm, Ytrain.hat_rf)
      Y1test.hats <- cbind(Y1test.hat_glm_1, Y1test.hat_glm_2, Y1test.hat_glm_3, Y1test.hat_dl, Y1test.hat_gbm, Y1test.hat_rf)
      Y0test.hats <- cbind(Y0test.hat_glm_1, Y0test.hat_glm_2, Y0test.hat_glm_3, Y0test.hat_dl, Y0test.hat_gbm, Y0test.hat_rf)

      # nnls
      fit.nnlsOR <- nnls(A=Ytrain.hats, b=Y)
      initCoefOR <- coef(fit.nnlsOR)
      initCoefOR[is.na(initCoefOR)] <- 0.0

      # normalize so sum(coef) = 1 if possible
      if (sum(initCoefOR) > 0) {
        coefOR <- initCoefOR / sum(initCoefOR)
      } else {
        warning("All algorithms have zero weight", call. = FALSE)
        coefOR <- initCoefOR
      }

      Y1test.comb <- as.numeric(crossprod(t(Y1test.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))
      Y0test.comb <- as.numeric(crossprod(t(Y0test.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))

      Yall[testIndex] = Y[testIndex]
      Y1hatall[testIndex] = Y1test.comb
      Y0hatall[testIndex] = Y0test.comb
      Zall[testIndex] = Z[testIndex]
      Zhatall[testIndex] = Ztest.comb
      dataALL[testIndex, ] = data[testIndex, ]
    }

    estimates_DML2 = DR(Y=Yall,Z=Zall,interZ=interZ, data=dataALL,
                        Z.hat=Zhatall,Y1.hat = Y1hatall, Y0.hat = Y0hatall)
    overallEst[m,] = as.numeric(estimates_DML2[,1])
    overallSE[m,] = as.numeric(estimates_DML2[,2])

    print(paste(m,"round out of ",mCrossFit," and estimated values:",
                paste(round(overallEst[m,],3),collapse=",")))
    print(paste(m,"round out of ",mCrossFit," and estimated SE:",
                paste(round(overallSE[m,],3),collapse=",")))

  }

  Estimate = apply(overallEst,2,median)
  SE = apply(overallSE,2,median)
  tau.est <- cbind(Estimate=Estimate,SE=SE)

  colnames(tau.est) <-  c("Estimate", "Std. Error")
  rownames(tau.est) <- colnames(interZ.mat)

  Z.hats = cbind(Z.hat=Ztest.comb, Ztest.hats)
  Y1.hats = cbind(Y1.hat=Y1test.comb, Y1test.hats)
  Y0.hats = cbind(Y0.hat=Y0test.comb, Y0test.hats)

  if (!crossFitting) {
    ans <- list(coef.ER=coefER, coef.OR=coefOR,
                Estimate=tau.est,
                Z.hat = Ztest.comb,  Y1.hat = Y1test.comb, Y0.hat = Y0test.comb,
                Z.hats = data.frame(Ztest.hats), Y1.hats = data.frame(Y1test.hats), Y0.hats = data.frame(Y0test.hats))
  } else {
    ans <- list(Estimate=tau.est)
  }

  class(ans) <- "CURobustML"

  return(ans)
}


#' Double demeaning estimator with treatment and outcome learners
#'
#' @param Y continuous outcome variable
#' @param Z binary treatment indicator, 1 - treatment, 0 - control
#' @param X vector, matrix, or dataframe containing measured confounders
#' @param interZ formula that contains the variables that "interact"  with the treatment. "1" will be always added. The default is no interaction, i.e., formula = formula(~1).
#' @param ID cluster identifier
#' @param data dataframe containing the variables in the model
#' @param library character vector of prediction algorithms. The available methods are  glm, deeplearning, gbm, and randomForests. The default methods are glm and deeplearning.
#' @param crossFitting whether to do cross fitting. The default is FALSE, and currently it is not available.
#' @param K number of folds. The default is 5, and currently it is not available.
#' @param mCrossFit number of cross fitting. The default is FALSE, and currently it is not available.
#'
#' @return
#' An \code{DDcomb} with the following elements:
#'    \item{coef.ER}{vector of the coefficients for prediction algorithms in the demeaned treatment model}
#'    \item{coef.OR}{vector of the coefficients for prediction algorithms in the demeaned outcome model}
#'    \item{Estimate}{estimates and standard errors of treatment effects}
#'    \item{Z.hat}{final weighted prediction for the demeaned treatment}
#'    \item{Y1.hat}{final weighted prediction for the demeaned outcome among treated units}
#'    \item{Y0.hat}{final weighted prediction for the demeaned outcome among control units}
#'    \item{Z.hats}{all the predictions for the demeaned treatment from prediction algorithms}
#'    \item{Y1.hats}{all the predictions for the demeaned outcome among treated units from prediction algorithms}
#'    \item{Y0.hats}{all the predictions for the demeaned outcome among control units from prediction algorithms}
#' @export
#'
#' @importFrom stats ave formula median model.matrix predict
#' @importFrom nnls nnls
#' @importFrom caret createFolds
#' @importFrom h2o as.h2o h2o.glm h2o.gbm h2o.deeplearning h2o.randomForest
#' @examples
#' # with two-level data
#' DDcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
#'  X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)
#'
#' \donttest{
#' # with cross-classified data
#' DDcomb(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1),
#'  X=crossclassified_data[, c("X1", "X2", "X3", "W1", "Q1")], ID=crossclassified_data$f12id,
#'  data=crossclassified_data)
#' }
DDcomb <- function(Y, Z, X, interZ=formula(~1), ID, data, library=c("glm", "deeplearning"), crossFitting=FALSE,
                   K = 5, mCrossFit = 100) {

  if(length(library) != sum(library %in% c("glm", "deeplearning", "gbm", "randomForests"))) {
    stop("library is not in the appropriate format. Choose methods among glm, gbm, randomForests, and deeplearning.")
  }

  if (crossFitting) {
    stop("Cross fitting has not yet been implemented.")
  }

  if(!crossFitting) {
    print("No crossFitting declared. Overriding mCrossfit and K variables to 1 both.")
    mCrossFit = 1; K = 1
  }

  # variables
  IDfactor = as.factor(ID)
  IDnum = as.numeric(IDfactor)

  # group deviation & group mean
  Y_demeaned = Y - ave(Y,IDnum )
  Z_demeaned = Z - ave(Z,IDnum )
  Z_mean = ave(Z,IDnum )
  X_demeaned = apply(as.matrix(X), 2, function(y) y - ave(y, IDnum))

  # Data structures for h2o
  lvl2var <-  colSums(apply(as.matrix(X), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1 <- as.matrix(X)[, !lvl2var]
  X_demeaned_lvl1 <-  as.matrix(X_demeaned)[, !lvl2var]
  X_demeaned_lvl1_num <- as.matrix(as.matrix(X_demeaned_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])

  X_lvl1_mean <- apply(as.matrix(X)[, !lvl2var], 2, function(y) ave(y, IDnum))
  W_lvl2 <- as.matrix(X)[, lvl2var]

  # datasets
  h2oDataE = as.h2o(data.frame(Z_demeaned,X_demeaned_lvl1, X_lvl1_mean, W_lvl2, X_demeaned_lvl1_num^2))

  h2oAll = as.h2o(data.frame(Y_demeaned,Z_demeaned,X_demeaned_lvl1, X_lvl1_mean, W_lvl2, X_demeaned_lvl1_num^2))
  h2oAll_Z1 = h2oAll_Z0 = h2oAll
  h2oAll_Z1$Z_demeaned = as.h2o(rep(1,length(Y_demeaned)) - ave(Z,IDnum)); h2oAll_Z0$Z_demeaned = as.h2o(rep(0,length(Y_demeaned)) - ave(Z,IDnum))

  interZ.mat <- model.matrix(interZ, data=data)

  if(is.null(interZ)) {
    overallEst = overallSE = matrix(0,mCrossFit,1)
  } else {
    overallEst = overallSE = matrix(0,mCrossFit,ncol(interZ.mat))
  }


  for(m in 1:mCrossFit) {

    dataALL = data.frame(matrix(0, nrow(data), ncol(data)))
    colnames(dataALL) = colnames(data)
    Yall = rep(0,length(Y))
    Y1hatall = rep(0,length(Y))
    Y0hatall = rep(0,length(Y))
    Zall = rep(0,length(Z))
    Zhatall = rep(0,length(Z))

    splitIndex = createFolds(IDfactor,k = K,list=FALSE)

    if(is.null(interZ)) {
      overallFold = matrix(0,K,1)
    } else {
      overallFold = matrix(0,K,ncol(interZ.mat))
    }


    for(i in 1:K) {
      trainIndex = which(splitIndex == i)
      if(!crossFitting) {
        testIndex = trainIndex
      } else {
        testIndex = which(splitIndex != i)
      }

      # ::: propensity score model
      h2oDataETrain = h2oDataE[trainIndex,]; h2oDataETest = h2oDataE[testIndex,]
      h2oAllTrain = h2oAll[trainIndex,]; h2oAllTest = h2oAll[testIndex,]
      h2oAllTest_Z1 = h2oAll_Z1[testIndex,]; h2oAllTest_Z0 = h2oAll_Z0[testIndex,]

      fitER_glm <- fitER_dl <- fitER_gbm <- fitER_rf <- Ztest.hat_glm <- Ztest.hat_dl <- Ztest.hat_gbm <- Ztest.hat_rf <-
        Ztrain.hat_glm <- Ztrain.hat_dl <- Ztrain.hat_gbm <- Ztrain.hat_rf <- NULL

      if ("glm" %in% library) {
        fitER_glm = h2o.glm(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataETrain,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oDataE) - ncol(X_demeaned_lvl1_num)))
        Ztrain.hat_glm = as.numeric(as.data.frame(predict(fitER_glm,h2oDataETrain)$pred)$predict)
        Ztest.hat_glm = as.numeric(as.data.frame(predict(fitER_glm,h2oDataETest)$pred)$predict)
      }

      if ("deeplearning" %in% library) {
        fitER_dl = h2o.deeplearning(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataETrain,hidden= rep((ncol(h2oDataE) -2) + nlevels(IDfactor),2),epochs=100,adaptive_rate=FALSE,rate=10^-5)
        Ztrain.hat_dl = as.numeric(as.data.frame(predict(fitER_dl,h2oDataETrain)$pred)$predict)
        Ztest.hat_dl = as.numeric(as.data.frame(predict(fitER_dl,h2oDataETest)$pred)$predict)
      }

      if ("gbm" %in% library) {
        fitER_gbm = h2o.gbm(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataETrain,ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
        Ztrain.hat_gbm = as.numeric(as.data.frame(predict(fitER_gbm,h2oDataETrain)$pred)$predict)
        Ztest.hat_gbm = as.numeric(as.data.frame(predict(fitER_gbm,h2oDataETest)$pred)$predict)
      }

      if ("randomForests" %in% library) {
        fitER_rf = h2o.randomForest(x=2:ncol(h2oDataE),y=1,training_frame=h2oDataETrain,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
        Ztrain.hat_rf = as.numeric(as.data.frame(predict(fitER_rf,h2oDataETrain)$pred)$predict)
        Ztest.hat_rf = as.numeric(as.data.frame(predict(fitER_rf,h2oDataETest)$pred)$predict)
      }

      # all predictions
      Ztrain.hats <- cbind(Ztrain.hat_glm, Ztrain.hat_dl, Ztrain.hat_gbm, Ztrain.hat_rf)
      Ztest.hats <- cbind(Ztest.hat_glm, Ztest.hat_dl, Ztest.hat_gbm, Ztest.hat_rf)

      # nnls
      fit.nnlsER <- nnls(A=Ztrain.hats, b=Z)
      initCoefER <- coef(fit.nnlsER)
      initCoefER[is.na(initCoefER)] <- 0.0

      # normalize so sum(coef) = 1 if possible
      if (sum(initCoefER) > 0) {
        coefER <- initCoefER / sum(initCoefER)
      } else {
        warning("All algorithms have zero weight", call. = FALSE)
        coefER <- initCoefER
      }

      Ztest.comb <- as.numeric(crossprod(t(Ztest.hats[, coefER != 0, drop = FALSE]), coefER[coefER != 0]))

      # ::: outcome model
      fitOR_glm <- fitOR_dl <- fitOR_gbm <- fitOR_rf <- Ytrain.hat_glm <- Ytrain.hat_dl <- Ytrain.hat_gbm <- Ytrain.hat_rf <-
        Y1test.hat_glm <- Y1test.hat_dl <- Y1test.hat_gbm <- Y1test.hat_rf <- Y0test.hat_glm <- Y0test.hat_dl <- Y0test.hat_gbm <- Y0test.hat_rf <- NULL

      if ("glm" %in% library) {
        fitOR_glm = h2o.glm(x=2:ncol(h2oAll),y=1,training_frame=h2oAllTrain,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAll) -ncol(X_demeaned_lvl1_num)))
        Ytrain.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAllTrain)$pred)$predict)
        Y1test.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_glm = as.numeric(as.data.frame(predict(fitOR_glm,h2oAllTest_Z0)$pred)$predict)
      }

      if ("deeplearning" %in% library) {
        fitOR_dl = h2o.deeplearning(x=2:ncol(h2oAll),y=1,training_frame=h2oAllTrain,hidden= rep((ncol(h2oAll) -2) + nlevels(IDfactor),2),epochs=100,adaptive_rate=FALSE,rate=10^-5)
        Ytrain.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTrain)$pred)$predict)
        Y1test.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTest_Z0)$pred)$predict)
      }

      if ("gbm" %in% library) {
        fitOR_gbm = h2o.gbm(x=2:ncol(h2oAll),y=1,training_frame=h2oAllTrain, ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
        Ytrain.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTrain)$pred)$predict)
        Y1test.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTest_Z0)$pred)$predict)
      }

      if ("randomForests" %in% library) {
        fitOR_rf = h2o.randomForest(x=2:ncol(h2oAll),y=1,training_frame=h2oAllTrain,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
        Ytrain.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTrain)$pred)$predict)
        Y1test.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTest_Z0)$pred)$predict)
      }

      # all predictions
      Ytrain.hats <- cbind(Ytrain.hat_glm, Ytrain.hat_dl, Ytrain.hat_gbm, Ytrain.hat_rf)
      Y1test.hats <- cbind(Y1test.hat_glm, Y1test.hat_dl, Y1test.hat_gbm, Y1test.hat_rf)
      Y0test.hats <- cbind(Y0test.hat_glm, Y0test.hat_dl, Y0test.hat_gbm, Y0test.hat_rf)

      # nnls
      fit.nnlsOR <- nnls(A=Ytrain.hats, b=Y)
      initCoefOR <- coef(fit.nnlsOR)
      initCoefOR[is.na(initCoefOR)] <- 0.0

      # normalize so sum(coef) = 1 if possible
      if (sum(initCoefOR) > 0) {
        coefOR <- initCoefOR / sum(initCoefOR)
      } else {
        warning("All algorithms have zero weight", call. = FALSE)
        coefOR <- initCoefOR
      }

      Y1test.comb <- as.numeric(crossprod(t(Y1test.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))
      Y0test.comb <- as.numeric(crossprod(t(Y0test.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))

      Yall[testIndex] = Y[testIndex]
      Y1hatall[testIndex] = Y1test.comb
      Y0hatall[testIndex] = Y0test.comb
      Zall[testIndex] = Z[testIndex]
      Zhatall[testIndex] = Ztest.comb
      dataALL[testIndex, ] = data[testIndex, ]
      IDnumtest <- IDnum[testIndex]
    }

    estimates_DML2 = DD(Y=Yall,Z=Zall,interZ=interZ, data=dataALL, ID=IDnumtest,
                        Z.hat=Zhatall, Y0.hat = Y0hatall)
    overallEst[m,] = as.numeric(estimates_DML2[,1])
    overallSE[m,] = as.numeric(estimates_DML2[,2])

    print(paste(m,"round out of ",mCrossFit," and estimated values:",
                paste(round(overallEst[m,],3),collapse=",")))
    print(paste(m,"round out of ",mCrossFit," and estimated SE:",
                paste(round(overallSE[m,],3),collapse=",")))

  }

  Estimate = apply(overallEst,2,median)
  SE = apply(overallSE,2,median)
  tau.est <- cbind(Estimate=Estimate,SE=SE)

  colnames(tau.est) <-  c("Estimate", "Std. Error")
  rownames(tau.est) <- colnames(interZ.mat)

  Z.hats = cbind(Z.hat=Ztest.comb, Z.hat_glm=Ztest.hat_glm, Z.hat_dl=Ztest.hat_dl, Z.hat_gbm=Ztest.hat_gbm, Z.hat_rf=Ztest.hat_rf)
  Y1.hats = cbind(Y1.hat=Y1test.comb, Y1.hat_glm=Y1test.hat_glm, Y1.hat_dl=Y1test.hat_dl, Y1.hat_gbm=Y1test.hat_gbm, Y1.hat_rf=Y1test.hat_rf)
  Y0.hats = cbind(Y0.hat=Y0test.comb, Y0.hat_glm=Y0test.hat_glm, Y0.hat_dl=Y0test.hat_dl, Y0.hat_gbm=Y0test.hat_gbm, Y0.hat_rf=Y0test.hat_rf)

  if (!crossFitting) {
    ans <- list(coef.ER=coefER, coef.OR=coefOR,
                Estimate=tau.est,
                Z.hat = Ztest.comb,  Y1.hat = Y1test.comb, Y0.hat = Y0test.comb,
                Z.hats = data.frame(Ztest.hats),
                Y1.hats = data.frame(Y1test.hats), Y0.hats = data.frame(Y0test.hats))
  } else {
    ans <- list(Estimate=tau.est)
  }

  class(ans) <- "CURobustML"

  return(ans)
}


#' Double demeaning estimator with proxy regressions with treatment and outcome learners
#'
#' @param Y continuous outcome variable
#' @param Z binary treatment indicator, 1 - treatment, 0 - control
#' @param X vector, matrix, or dataframe containing measured confounders
#' @param interZ formula that contains the variables that "interact"  with the treatment. "1" will be always added. The default is no interaction, i.e., formula = formula(~1).
#' @param ID cluster identifier
#' @param data dataframe containing the variables in the model
#' @param library character vector of prediction algorithms. The available methods are  glm, deeplearning, gbm, and randomForests. The default methods are glm and deeplearning. Three types of glm include: glm with cluster dummies, glm with cluster-constant components of individual-level covariates, and glmm with random effects of clusters.
#' @param crossFitting whether to do cross fitting. The default is FALSE, and currently it is not available.
#' @param K number of folds. The default is 5, and currently it is not available.
#' @param mCrossFit number of cross fitting. The default is FALSE, and currently it is not available.
#'
#' @return
#' An \code{DDPRcomb} with the following elements:
#'    \item{coef.ER}{vector of the coefficients for prediction algorithms in the treatment model}
#'    \item{coef.OR}{vector of the coefficients for prediction algorithms in the outcome model}
#'    \item{Estimate}{estimates and standard errors of treatment effects}
#'    \item{Z.hat}{final weighted prediction for the treatment}
#'    \item{Y1.hat}{final weighted prediction for the outcome among treated units}
#'    \item{Y0.hat}{final weighted prediction for the outcome among control units}
#'    \item{Z.hats}{all the predictions for the treatment from prediction algorithms}
#'    \item{Y1.hats}{all the predictions for the outcome among treated units from prediction algorithms}
#'    \item{Y0.hats}{all the predictions for the outcome among control units from prediction algorithms}
#'
#' @export
#'
#' @importFrom stats ave formula median model.matrix predict
#' @importFrom nnls nnls
#' @importFrom caret createFolds
#' @importFrom WeightIt make_full_rank
#' @importFrom lme4 glmer lmer glmerControl lmerControl
#' @importFrom h2o as.h2o h2o.glm h2o.gbm h2o.deeplearning h2o.randomForest
#' @examples
#' # with two-level data
#' DDPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
#'  X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)
#'
#' \donttest{
#' # with cross-classified data
#' DDPRcomb(Y=crossclassified_data$Y, Z=crossclassified_data$Z, interZ=(~ W1),
#'  X=crossclassified_data[, c("X1", "X2", "X3", "W1", "Q1")], ID=crossclassified_data$f12id,
#'  data=crossclassified_data)
#' }
DDPRcomb <- function(Y, Z, X, interZ=formula(~1), ID, data, library=c("glm", "deeplearning"), crossFitting=FALSE,
                     K = 5, mCrossFit = 100) {

  if(length(library) != sum(library %in% c("glm", "deeplearning", "gbm", "randomForests"))) {
    stop("library is not in the appropriate format. Choose methods among glm, gbm, randomForests, and deeplearning.")
  }

  if (crossFitting) {
    stop("Cross fitting has not yet been implemented.")
  }

  if(!crossFitting) {
    print("No crossFitting declared. Overriding mCrossfit and K variables to 1 both.")
    mCrossFit = 1; K = 1
  }

  # variables
  IDfactor = as.factor(ID)
  IDnum = as.numeric(IDfactor)
  Zfactor = as.factor(Z)
  Z = as.numeric(as.character(Z))

  lvl2var <- colSums(apply(as.matrix(X), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1 <- as.matrix(X)[, !lvl2var]
  W_lvl2 <- as.matrix(X)[, lvl2var]
  X_lvl1_num <- as.matrix(as.matrix(X_lvl1)[, apply(as.matrix(X_lvl1), 2, function(x) sum(!unique(x) %in% c(0, 1))) > 0])
  X_lvl1_mean <- apply(as.matrix(X)[, !lvl2var], 2, function(y) ave(y, IDnum))

  X_lvl1_num_sq <- X_lvl1_num^2
  colnames(X_lvl1_num_sq) <- paste0(colnames(X_lvl1_num_sq), "_sq")

  # datasets for the propensity scores
  h2oDataE= as.h2o(data.frame(Zfactor, X_lvl1, W_lvl2, IDfactor, X_lvl1_num_sq))
  h2oDataE2 = as.h2o(data.frame(Zfactor, X_lvl1, W_lvl2, X_lvl1_mean, X_lvl1_num_sq))

  VarAll = data.frame(X_lvl1,W_lvl2)
  VarInt.m <- model.matrix(~ .^2, data=VarAll)
  VarInt.m2 <- WeightIt::make_full_rank(VarInt.m, with.intercept = TRUE)
  REdataE = data.frame(VarInt.m2, X_lvl1_num_sq, Zfactor, IDfactor)

  covFrml <- paste(colnames(VarInt.m2), collapse="+")
  glmerFrml <- formula(paste("Zfactor ~", covFrml, " + (1 | IDfactor)"))

  # datasets for the outcome
  h2oAll = as.h2o(data.frame(Y, Z, X_lvl1, W_lvl2 ,IDfactor, X_lvl1_num_sq))
  h2oAll_Z1 = h2oAll_Z0 = h2oAll; h2oAll_Z1$Z = 1; h2oAll_Z0$Z = 0

  h2oAll2 = as.h2o(data.frame(Y, Z, X_lvl1, W_lvl2, X_lvl1_mean, X_lvl1_num_sq))
  h2oAll2_Z1 = h2oAll2_Z0 = h2oAll2; h2oAll2_Z1$Z = 1; h2oAll2_Z0$Z = 0

  VarTmp = data.frame(Z,X_lvl1,W_lvl2)
  VarAll2 = rbind(VarTmp,VarTmp,VarTmp)
  VarAll2$Z[(nrow(VarTmp) + 1):(2*nrow(VarTmp))] = 1
  VarAll2$Z[(2*nrow(VarTmp) + 1):(3*nrow(VarTmp))] = 0
  VarInt2.m <- model.matrix(~ .^2, data=VarAll2)
  VarInt2.m2 <- WeightIt::make_full_rank(VarInt2.m, with.intercept = TRUE)

  Var2Tmp = data.frame(Y, X_lvl1_num_sq, IDfactor)
  REdata = data.frame(VarInt2.m2[1:nrow(VarTmp), ], Var2Tmp)
  REdata_Z1 = data.frame(VarInt2.m2[(nrow(VarTmp) + 1):(2*nrow(VarTmp)), ], Var2Tmp)
  REdata_Z0 = data.frame(VarInt2.m2[(2*nrow(VarTmp) + 1):(3*nrow(VarTmp)), ], Var2Tmp)

  covFrml2 <- paste(colnames(VarInt2.m2), collapse="+")
  lmerFrml <- formula(paste("Y ~", covFrml2, " + (1 | IDfactor)"))

  interZ.mat <- model.matrix(interZ, data=data)

  if(is.null(interZ)) {
    overallEst = overallSE = matrix(0,mCrossFit,1)
  } else {
    overallEst = overallSE = matrix(0,mCrossFit,ncol(interZ.mat))
  }


  for(m in 1:mCrossFit) {

    dataALL = data.frame(matrix(0, nrow(data), ncol(data)))
    colnames(dataALL) = colnames(data)
    Yall = rep(0,length(Y))
    Y1hatall = rep(0,length(Y))
    Y0hatall = rep(0,length(Y))
    Zall = rep(0,length(Z))
    Zhatall = rep(0,length(Z))

    splitIndex = createFolds(IDfactor,k = K,list=FALSE)

    if(is.null(interZ)) {
      overallFold = matrix(0,K,1)
    } else {
      overallFold = matrix(0,K,ncol(interZ.mat))
    }


    for(i in 1:K) {
      trainIndex = which(splitIndex == i)
      if(!crossFitting) {
        testIndex = trainIndex
      } else {
        testIndex = which(splitIndex != i)
      }

      # ::: propensity score model
      h2oDataETrain = h2oDataE[trainIndex,]; h2oDataETest = h2oDataE[testIndex,]
      h2oDataE2Train = h2oDataE2[trainIndex,]; h2oDataE2Test = h2oDataE2[testIndex,]

      h2oAllTrain = h2oAll[trainIndex,]; h2oAllTest_Z1 = h2oAll_Z1[testIndex,]; h2oAllTest_Z0 = h2oAll_Z0[testIndex,]
      h2oAll2Train = h2oAll2[trainIndex,]; h2oAll2Test_Z1 = h2oAll2_Z1[testIndex,]; h2oAll2Test_Z0 = h2oAll2_Z0[testIndex,]

      REdataETrain = REdataE[trainIndex,]; REdataETest = REdataE[testIndex,]
      REdataTrain = REdata[trainIndex,]; REdataTest_Z1 = REdata_Z1[testIndex,]; REdataTest_Z0 = REdata_Z0[testIndex,]

      fitER_glm <- fitER_dl <- fitER_gbm <- fitER_rf <- Ztest.hat_glm <- Ztest.hat_dl <- Ztest.hat_gbm <- Ztest.hat_rf <-
        Ztrain.hat_glm <- Ztrain.hat_dl <- Ztrain.hat_gbm <- Ztrain.hat_rf <- NULL

      if ("glm" %in% library) {
        # ensemble glm for propensity score training
        if(length(2:(ncol(h2oDataETrain) -ncol(X_lvl1_num)-1))>=2) {
          fitER_glm_1 = h2o.glm(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,alpha=1,lambda=0,family="binomial",interactions=2:(ncol(h2oDataETrain) -ncol(X_lvl1_num)-1))
          fitER_glm_2 = h2o.glm(x=2:ncol(h2oDataE2Train),y=1,training_frame=h2oDataE2Train,alpha=1,lambda=0,family="binomial",interactions=2:(ncol(h2oDataE2Train) -ncol(X_lvl1_num)))
        } else {
          fitER_glm_1 = h2o.glm(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,alpha=1,lambda=0,family="binomial")
          fitER_glm_2 = h2o.glm(x=2:ncol(h2oDataE2Train),y=1,training_frame=h2oDataE2Train,alpha=1,lambda=0,family="binomial")
        }

        fitER_glm_3 = glmer(glmerFrml, family=binomial(), data=REdataETrain, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5), calc.derivs = FALSE), nAGQ = 0)

        Ztrain.hat_glm_1 = as.numeric(as.data.frame(predict(fitER_glm_1,h2oDataETrain)[,"p1"])$p1)
        Ztrain.hat_glm_2 = as.numeric(as.data.frame(predict(fitER_glm_2,h2oDataE2Train)[,"p1"])$p1)
        Ztrain.hat_glm_3 = as.numeric(predict(fitER_glm_3, type="response"))

        Ztest.hat_glm_1 = as.numeric(as.data.frame(predict(fitER_glm_1,h2oDataETest)[,"p1"])$p1)
        Ztest.hat_glm_2 = as.numeric(as.data.frame(predict(fitER_glm_2,h2oDataE2Test)[,"p1"])$p1)
        Ztest.hat_glm_3 = as.numeric(predict(fitER_glm_3, type="response", REdataETest))

      }

      if ("deeplearning" %in% library) {
        fitER_dl = h2o.deeplearning(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,hidden= rep((ncol(h2oDataETrain) -2) + nlevels(IDfactor),2),epochs=100,adaptive_rate=FALSE,rate=10^-5)
        Ztrain.hat_dl = as.numeric(as.data.frame(predict(fitER_dl,h2oDataETrain)[,"p1"])$p1)
        Ztest.hat_dl = as.numeric(as.data.frame(predict(fitER_dl,h2oDataETest)[,"p1"])$p1)
      }

      if ("gbm" %in% library) {
        fitER_gbm = h2o.gbm(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
        Ztrain.hat_gbm = as.numeric(as.data.frame(predict(fitER_gbm,h2oDataETrain)$pred)$predict)
        Ztest.hat_gbm = as.numeric(as.data.frame(predict(fitER_gbm,h2oDataETest)$pred)$predict)
      }

      if ("randomForests" %in% library) {
        fitER_rf = h2o.randomForest(x=2:ncol(h2oDataETrain),y=1,training_frame=h2oDataETrain,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
        Ztrain.hat_rf = as.numeric(as.data.frame(predict(fitER_rf,h2oDataETrain)$pred)$predict)
        Ztest.hat_rf = as.numeric(as.data.frame(predict(fitER_rf,h2oDataETrain)$pred)$predict)
      }


      # all predictions
      Ztrain.hats <- cbind(Ztrain.hat_glm_1, Ztrain.hat_glm_2, Ztrain.hat_glm_3, Ztrain.hat_dl, Ztrain.hat_gbm, Ztrain.hat_rf)
      Ztest.hats <- cbind(Ztest.hat_glm_1, Ztest.hat_glm_2, Ztest.hat_glm_3, Ztest.hat_dl, Ztest.hat_gbm, Ztest.hat_rf)

      # nnls
      fit.nnlsER <- nnls(A=Ztrain.hats, b=Z)
      initCoefER <- coef(fit.nnlsER)
      initCoefER[is.na(initCoefER)] <- 0.0

      # normalize so sum(coef) = 1 if possible
      if (sum(initCoefER) > 0) {
        coefER <- initCoefER / sum(initCoefER)
      } else {
        warning("All algorithms have zero weight", call. = FALSE)
        coefER <- initCoefER
      }

      Ztest.comb <- as.numeric(crossprod(t(Ztest.hats[, coefER != 0, drop = FALSE]), coefER[coefER != 0]))

      # ::: outcome model
      fitOR_glm <- fitOR_dl <- fitOR_gbm <- fitOR_rf <- Ytrain.hat_glm <- Ytrain.hat_dl <- Ytrain.hat_gbm <- Ytrain.hat_rf <-
        Y1test.hat_glm <- Y1test.hat_dl <- Y1test.hat_gbm <- Y1test.hat_rf <- Y0test.hat_glm <- Y0test.hat_dl <- Y0test.hat_gbm <- Y0test.hat_rf <- NULL

      if ("glm" %in% library) {
        # ensemble glm for outcome reg training
        if(length(2:(ncol(h2oAllTrain) -ncol(X_lvl1_num)-1))>=2) {
          fitOR_glm_1 = h2o.glm(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAllTrain) -ncol(X_lvl1_num)-1))
          fitOR_glm_2 = h2o.glm(x=2:ncol(h2oAll2Train),y=1,training_frame=h2oAll2Train,alpha=1,lambda=0,family="gaussian",interactions=2:(ncol(h2oAll2Train) -ncol(X_lvl1_num)))
        } else {
          fitOR_glm_1 = h2o.glm(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain,alpha=1,lambda=0,family="gaussian")
          fitOR_glm_2 = h2o.glm(x=2:ncol(h2oAll2Train),y=1,training_frame=h2oAll2Train,alpha=1,lambda=0,family="gaussian")
        }

        fitOR_glm_3 = lmer(lmerFrml, data=REdataTrain, control=lmerControl(calc.derivs = FALSE))

        Ytrain.hat_glm_1 = as.numeric(as.data.frame(predict(fitOR_glm_1,h2oAllTrain)$pred)$predict)
        Ytrain.hat_glm_2 = as.numeric(as.data.frame(predict(fitOR_glm_2,h2oAll2Train)$pred)$predict)
        Ytrain.hat_glm_3 = as.numeric(predict(fitOR_glm_3))

        Y1test.hat_glm_1 = as.numeric(as.data.frame(predict(fitOR_glm_1,h2oAllTest_Z1)$pred)$predict)
        Y1test.hat_glm_2 = as.numeric(as.data.frame(predict(fitOR_glm_2,h2oAll2Test_Z1)$pred)$predict)
        Y1test.hat_glm_3 = as.numeric(predict(fitOR_glm_3, REdataTest_Z1))

        Y0test.hat_glm_1 = as.numeric(as.data.frame(predict(fitOR_glm_1,h2oAllTest_Z0)$pred)$predict)
        Y0test.hat_glm_2 = as.numeric(as.data.frame(predict(fitOR_glm_2,h2oAll2Test_Z0)$pred)$predict)
        Y0test.hat_glm_3 = as.numeric(predict(fitOR_glm_3, REdataTest_Z0))

      }

      if ("deeplearning" %in% library) {
        fitOR_dl = h2o.deeplearning(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain,hidden= rep((ncol(h2oAllTrain) -2) + nlevels(IDfactor),2), epochs=100,adaptive_rate=FALSE,rate=10^-5)
        Ytrain.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTrain)$pred)$predict)
        Y1test.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_dl = as.numeric(as.data.frame(predict(fitOR_dl,h2oAllTest_Z0)$pred)$predict)
      }

      if ("gbm" %in% library) {
        fitOR_gbm = h2o.gbm(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain, ntrees=8,max_depth=5,min_rows=min(5,min(table(IDfactor))/2))
        Ytrain.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTrain)$pred)$predict)
        Y1test.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_gbm = as.numeric(as.data.frame(predict(fitOR_gbm,h2oAllTest_Z0)$pred)$predict)
      }

      if ("randomForests" %in% library) {
        fitOR_rf = h2o.randomForest(x=2:ncol(h2oAllTrain),y=1,training_frame=h2oAllTrain,ntrees = 500, max_depth = 50, stopping_rounds = 2, stopping_tolerance = 1e-2, score_each_iteration = T)
        Ytrain.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTrain)$pred)$predict)
        Y1test.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTest_Z1)$pred)$predict)
        Y0test.hat_rf = as.numeric(as.data.frame(predict(fitOR_rf,h2oAllTest_Z0)$pred)$predict)
      }

      # all predictions
      Ytrain.hats <- cbind(Ytrain.hat_glm_1, Ytrain.hat_glm_2, Ytrain.hat_glm_3, Ytrain.hat_dl, Ytrain.hat_gbm, Ytrain.hat_rf)
      Y1test.hats <- cbind(Y1test.hat_glm_1, Y1test.hat_glm_2, Y1test.hat_glm_3, Y1test.hat_dl, Y1test.hat_gbm, Y1test.hat_rf)
      Y0test.hats <- cbind(Y0test.hat_glm_1, Y0test.hat_glm_2, Y0test.hat_glm_3, Y0test.hat_dl, Y0test.hat_gbm, Y0test.hat_rf)

      # nnls
      fit.nnlsOR <- nnls(A=Ytrain.hats, b=Y)
      initCoefOR <- coef(fit.nnlsOR)
      initCoefOR[is.na(initCoefOR)] <- 0.0

      # normalize so sum(coef) = 1 if possible
      if (sum(initCoefOR) > 0) {
        coefOR <- initCoefOR / sum(initCoefOR)
      } else {
        warning("All algorithms have zero weight", call. = FALSE)
        coefOR <- initCoefOR
      }

      Y1test.comb <- as.numeric(crossprod(t(Y1test.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))
      Y0test.comb <- as.numeric(crossprod(t(Y0test.hats[, coefOR != 0, drop = FALSE]), coefOR[coefOR != 0]))

      Yall[testIndex] = Y[testIndex]
      Y1hatall[testIndex] = Y1test.comb
      Y0hatall[testIndex] = Y0test.comb
      Zall[testIndex] = Z[testIndex]
      Zhatall[testIndex] = Ztest.comb
      dataALL[testIndex, ] = data[testIndex, ]
      IDnumtest <- IDnum[testIndex]
    }

    estimates_DML2 = DD(Y=Yall,Z=Zall,interZ=interZ, data=dataALL,
                        Z.hat=Zhatall, Y0.hat = Y0hatall, ID=IDnumtest)
    overallEst[m,] = as.numeric(estimates_DML2[,1])
    overallSE[m,] = as.numeric(estimates_DML2[,2])

    print(paste(m,"round out of ",mCrossFit," and estimated values:",
                paste(round(overallEst[m,],3),collapse=",")))
    print(paste(m,"round out of ",mCrossFit," and estimated SE:",
                paste(round(overallSE[m,],3),collapse=",")))

  }

  Estimate = apply(overallEst,2,median)
  SE = apply(overallSE,2,median)
  tau.est <- cbind(Estimate=Estimate,SE=SE)

  colnames(tau.est) <-  c("Estimate", "Std. Error")
  rownames(tau.est) <- colnames(interZ.mat)

  Z.hats = cbind(Z.hat=Ztest.comb, Ztest.hats)
  Y1.hats = cbind(Y1.hat=Y1test.comb, Y1test.hats)
  Y0.hats = cbind(Y0.hat=Y0test.comb, Y0test.hats)

  if (!crossFitting) {
    ans <- list(coef.ER=coefER, coef.OR=coefOR,
                Estimate=tau.est,
                Z.hat = Ztest.comb,  Y1.hat = Y1test.comb, Y0.hat = Y0test.comb,
                Z.hats = data.frame(Ztest.hats), Y1.hats = data.frame(Y1test.hats), Y0.hats = data.frame(Y0test.hats))
  } else {
    ans <- list(Estimate=tau.est)
  }

  class(ans) <- "CURobustML"

  return(ans)
}



#' Doubly robust estimator
#'
#' @param Y continuous outcome variable
#' @param Z binary treatment indicator, 1 - treatment, 0 - control
#' @param interZ formula that contains the variables that "interact"  with the treatment. "1" will be always added. The default is no interaction, i.e., formula = formula(~1).
#' @param data dataframe containing the variables in the model
#' @param Z.hat treatment/propensity score prediction from \code{\link{DRPRcomb}}
#' @param Y1.hat outcome prediction among treated units from \code{\link{DRPRcomb}}
#' @param Y0.hat outcome prediction among untreated units from \code{\link{DRPRcomb}}
#'
#' @return estimates and standard errors
#' @export
#'
#' @importFrom stats  model.matrix coef vcov lm
#' @examples
#' DRPRcomb.rslt <- DRPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
#'  X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)
#'
#' # with final predictions
#' DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRPRcomb.rslt$Z.hat,
#'  Y1.hat=DRPRcomb.rslt$Y1.hat, Y0.hat=DRPRcomb.rslt$Y0.hat, data=twolevel_data)
#'
#' # with predictions from glm with fixed effects of clusters
#' DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRPRcomb.rslt$Z.hats$Ztest.hat_glm_1,
#'  Y1.hat=DRPRcomb.rslt$Y1.hats$Y1test.hat_glm_1, Y0.hat=DRPRcomb.rslt$Y0.hats$Y0test.hat_glm_1,
#'  data=twolevel_data)
#'
#' # with predictions from glm wtih random effects of clusters
#  DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRPRcomb.rslt$Z.hats$Ztest.hat_glm_3,
#   Y1.hat=DRPRcomb.rslt$Y1.hats$Y1test.hat_glm_3, Y0.hat=DRPRcomb.rslt$Y0.hats$Y0test.hat_glm_3,
#   data=twolevel_data)
DR <- function(Y, Z, interZ=formula(~1), Z.hat, Y1.hat, Y0.hat, data) {

  DRest.ind = (Z * (Y- Y1.hat) / Z.hat + Y1.hat) - ((1 -Z) * (Y - Y0.hat) / (1 - Z.hat) + Y0.hat)

  interZ.mat <- model.matrix(interZ, data=data)

  ATE_interZ_model = lm(DRest.ind ~ interZ.mat -1)
  SE_interZ = sqrt(diag(vcov(ATE_interZ_model)))
  ATE_interZ =coef(ATE_interZ_model)
  tau.est <- cbind(Estimate=ATE_interZ,SE=SE_interZ)

  colnames(tau.est) <-  c("Estimate", "Std. Error")
  rownames(tau.est) <- colnames(interZ.mat)

  return(tau.est)
}


#' Double demeaning estimator
#'
#' @param Y continuous outcome variable
#' @param Z binary treatment indicator, 1 - treatment, 0 - control
#' @param interZ formula that contains the variables that "interact"  with the treatment. "1" will be always added. The default is no interaction, i.e., formula = formula(~1).
#' @param ID cluster identifier
#' @param data dataframe containing the variables in the model
#' @param Z.hat treatment prediction from \code{\link{DDcomb}}, \code{\link{DRPRcomb}}, or \code{\link{DDPRcomb}}
#' @param Y0.hat outcome prediction among untreated units from \code{\link{DDcomb}}, \code{\link{DRPRcomb}}, or \code{\link{DDPRcomb}}
#'
#' @return estimates and standard errors
#' @export
#'
#' @importFrom stats ave model.matrix
#' @importFrom AER ivreg
#' @examples
#' DDcomb.rslt <- DDcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
#'  X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)
#'
#' # with final predictions
#' DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), ID=twolevel_data$id,
#'  Z.hat=DDcomb.rslt$Z.hat, Y0.hat=DDcomb.rslt$Y0.hat, data=twolevel_data)
#' # with predictions from glm
#' DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), ID=twolevel_data$id,
#'  Z.hat=DDcomb.rslt$Z.hats$Ztest.hat_glm, Y0.hat=DDcomb.rslt$Y0.hats$Y0test.hat_glm,
#'  data=twolevel_data)
#' # with predictions from deep learning
#' DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), ID=twolevel_data$id,
#'  Z.hat=DDcomb.rslt$Z.hats$Ztest.hat_dl, Y0.hat=DDcomb.rslt$Y0.hats$Y0test.hat_dl,
#'  data=twolevel_data)
#'
#' \donttest{
#' DDPRcomb.rslt <- DDPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
#'  X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)
#'
#' # with final predictions
#' DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), ID=twolevel_data$id,
#'  Z.hat=DDPRcomb.rslt$Z.hat, Y0.hat=DDPRcomb.rslt$Y0.hat, data=twolevel_data)
#' # with predictions from glm
#' DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), ID=twolevel_data$id,
#'  Z.hat=DDPRcomb.rslt$Z.hats$Ztest.hat_glm_1, Y0.hat=DDPRcomb.rslt$Y0.hats$Y0test.hat_glm_1,
#'  data=twolevel_data)
#' # with predictions from deep learning
#' DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), ID=twolevel_data$id,
#'  Z.hat=DDPRcomb.rslt$Z.hats$Ztest.hat_dl, Y0.hat=DDPRcomb.rslt$Y0.hats$Y0test.hat_dl,
#'  data=twolevel_data)
#' }
DD <- function(Y, Z, interZ=formula(~1),ID, Z.hat, Y0.hat, data) {

  Ygrpcent = Y - ave(Y, ID)
  Zgrpcent = Z - ave(Z, ID)

  Y.resid.D = Ygrpcent - Y0.hat
  Z.resid.D = Zgrpcent - Z.hat

  Y.resid.DD <- Y.resid.D - ave(Y.resid.D, ID)
  Z.resid.DD <- Z.resid.D - ave(Z.resid.D, ID)

  Zmat <- model.matrix(interZ, data=data)
  Zmat.rr <- apply(Zmat, 2, '*', Z.resid.DD)
  Zmat.Z <- apply(Zmat, 2, '*', Z)
  Zmat.Z.grpcent <- apply(Zmat.Z, 2, function(x) x -ave(x, ID))

  ivm = summary(ivreg(Y.resid.DD ~ Zmat.Z.grpcent -1 | Zmat.rr - 1))
  tau.est = matrix(0,ncol(Zmat),2)
  tau.est[,1] = ivm$coefficients[,1]
  tau.est[,2] = ivm$coefficients[,2]

  rownames(tau.est) = colnames(Zmat)
  colnames(tau.est) = c("Estimate","Std. Error")
  return(tau.est)
}



#' Summarize the estimates from the \code{CURobustML} routine
#'
#' @description coefficient table containing estimates, standard errors, z-value, and p-value. If exists, bootstrap standard errors are used by default.
#'
#' @param object object of class \code{CURobustML}.
#' @param digits integer indicating the number of decimal places
#' @param ... additional arguments ...
#'
#' @importFrom stats pnorm
#' @export
#'
summary.CURobustML <-
  function(object, digits = max(3L, getOption("digits") - 3L), ...){

    tau.est <- object$Estimate

    zvalue <- tau.est[,1] / tau.est[,2]
    pvalue <- 2 * pnorm(-abs(zvalue))

    coef.table <- as.matrix(cbind(tau.est[,1], tau.est[,2], zvalue, pvalue))

    colnames(coef.table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    rownames(coef.table) <- names(tau.est[,1])

    class(coef.table) <- "summary.CURobustML"

    return(coef.table)
  }



#' Print a coefficient table
#'
#' @param x object of class summary
#' @param digits integer indicating the number of decimal places
#' @param signif.stars whether to add significant starts
#' @param ... additional arguments ...
#'
#' @importFrom stats printCoefmat
#' @export
#'
print.summary.CURobustML <-
  function(x, digits = max(3L, getOption("digits") - 3L),
           signif.stars = getOption("show.signif.stars"), ...){
    cat("Coefficients: \n \n")

    printCoefmat(x, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)

  }


#' @title Setup function
#'
#' @description setup function to implement CURobustML
#'
#' @param ... additional arguments ...
#' @export
#'
#' @importFrom h2o h2o.init h2o.removeAll h2o.no_progress
#' @examples
#' CURobustML.setup()
CURobustML.setup <- function(...) {

  localH2O = h2o.init(...)
  h2o.removeAll()
  h2o.no_progress()

}

#' @title Overlap plot
#'
#' @description overlap plot between the treated group and untreated group
#'
#' @param X numeric vector (e.g., propensity score logit)
#' @param Z binary Treatment indicator, 1 - treatment, 0 - control
#' @param bin number of bins for histogram
#' @param ... additional arguments ...
#'
#' @return overlap plot between treated units and control units
#' @export
#'
#' @importFrom graphics axis hist par plot text
#' @examples
#' overlap(X=twolevel_data$lps, Z=twolevel_data$Z)
overlap <- function(X, Z, bin = 20, ...) {
  r1 <- range(X)
  if (!is.numeric(Z)) Z <- as.numeric(Z) - 1
  c.dat <- hist(X[Z == 0], seq(r1[1], r1[2], length = bin), plot = F)  # histogram data for control group
  t.dat <- hist(X[Z == 1], seq(r1[1], r1[2], length = bin), plot = F)  # histogram data for treatm. group
  c.dat$counts <- -c.dat$counts
  plot(t.dat, axes = F, ylim = c(min(c.dat$counts), max(t.dat$counts)), ...)
  plot(c.dat, add = T, density = 30)
  axis(1)
  ax.loc <- axis(2, labels = F)
  axis(2, at = ax.loc, labels = abs(ax.loc))
  y <- par('usr')[3:4]
  text(rep(max(X), 2), c(y[2] - diff(y)*.05, y[1] + diff(y)*.05),
       c('Treatment', 'Control'), adj = 1, xpd = T)
}


