context("testing main functions in CURobustML")

test_that("testing functions with glm only", {

  CURobustML.setup()

  DRPRcomb.rslt <- DRPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                            X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)

  DRPRrlst <- DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRPRcomb.rslt$Z.hats$Ztest.hat_glm_1,
                 Y1.hat=DRPRcomb.rslt$Y1.hats$Y1test.hat_glm_1, Y0.hat=DRPRcomb.rslt$Y0.hats$Y0test.hat_glm_1,
                 data=twolevel_data)

  expect_equal(DRPRrlst[,1], c(`(Intercept)` = 1.99338981595332, W1 = 0.957614872173895))

  DDcomb.rslt <- DDcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                        X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)

  DDrslt <- DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), ID=twolevel_data$id,
               Z.hat=DDcomb.rslt$Z.hats$Ztest.hat_glm, Y0.hat=DDcomb.rslt$Y0.hats$Y0test.hat_glm,
               data=twolevel_data)

  expect_equal(DDrslt[,1], c(`(Intercept)` = 2.03808160391043, W1 = 0.970151985228184))

  DDPRcomb.rslt <- DDPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                            X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)

  DDPRrslt <- DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), ID=twolevel_data$id,
                 Z.hat=DDPRcomb.rslt$Z.hats$Ztest.hat_glm_1, Y0.hat=DDPRcomb.rslt$Y0.hats$Y0test.hat_glm_1,
                 data=twolevel_data)

  expect_equal(DDPRrslt[,1], c(`(Intercept)` = 2.04076267488469, W1 = 0.957571008787558))

})
