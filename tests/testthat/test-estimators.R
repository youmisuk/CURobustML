context("Implementingestimators")

test_that("Whether DRPR gives us the same results", {

  CURobustML.setup()

  DRcomb.rslt <- DRPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                          X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)

  # with predictions from the first glm algorithm
  DRest <- DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DRcomb.rslt$Z.hats_glm$Ztest.hat_glm_1,
     Y1.hat=DRcomb.rslt$Y1.hats_glm$Y1test.hat_glm_1, Y0.hat=DRcomb.rslt$Y0.hats_glm$Y0test.hat_glm_1, data=twolevel_data)

  expect_equal(dput(DRest[1,1]), 1.99338981595332)

})

test_that("Whether DD gives us the same results", {

  CURobustML.setup()

  DDcomb.rslt <- DDcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                          X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)

  # with predictions from the first glm algorithm
  DDest <- DD(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DDcomb.rslt$Z.hats$Ztest.hat_glm,
              Y0.hat=DDcomb.rslt$Y0.hats$Y0test.hat_glm, data=twolevel_data, ID=twolevel_data$id)

  expect_equal(dput(DDest[1,1]), 2.03808160391043)

})

test_that("Whether DDPR gives us the same results", {

  CURobustML.setup()

  DDPRcomb.rslt <- DDPRcomb(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1),
                          X=twolevel_data[, c("X1", "X2", "X3", "W1")], ID=twolevel_data$id, data=twolevel_data)

  # with predictions from the first glm algorithm
  DDPRest <- DR(Y=twolevel_data$Y, Z=twolevel_data$Z, interZ=(~ W1), Z.hat=DDPRcomb.rslt$Z.hats_glm$Ztest.hat_glm_1,
              Y1.hat=DDPRcomb.rslt$Y1.hats_glm$Y1test.hat_glm_1, Y0.hat=DDPRcomb.rslt$Y0.hats_glm$Y0test.hat_glm_1, data=twolevel_data)

  expect_equal(dput(DDPRest[1,1]), 1.99338981595332)

})
