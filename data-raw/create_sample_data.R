# create sample data

set.seed(519)

twolevel.pop <- function(ids=1, Smpl.size = "150.30", trt.val=2, int.trt.val=1, sel.coef=0.3, out.coef=2,
                         intX2W1=FALSE, intX1X2=FALSE, intZW1=FALSE, intZX2=FALSE, intZX3=FALSE, intZX2sq=FALSE, intX1sq=FALSE) {

  # ::::: 1) generate the number of observations in each cluster :::::

  size <- strsplit(Smpl.size, "\\.")[[1]]  # the number of cluster and cluster sizes

  J <- as.numeric(size[1]) # level-2 unit, the number of cluster
  n.clus <- as.numeric(size[2]) # level-1 unit, cluster sizes

  if (n.clus < 20) {  # variation of cluster sizes
    sd.clus <- 0.5
  } else if (n.clus >= 20 & n.clus < 30) {
    sd.clus <- 1
  } else if (n.clus >= 30) {
    sd.clus <- 2
  }

  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes
  id <- as.factor(rep(ids:(J + ids - 1), N))              # cluster id

  # ::::: 2) generate a level-2 covariate, W :::::

  W1 <- rnorm(J, 0, 1)
  W2 <- rnorm(J, 0, 1)
  names(W1) <- names(W2)  <- levels(id)

  # ::::: 3) generate level-1 covariates, Xs with cluster-specific means :::::

  totalN <- length(id)
  X1 <- runif(totalN, -1, 1)
  X2 <- rnorm(totalN)
  X3_U <- runif(totalN, 0, 1)
  X3 <- ifelse(X3_U < 0.5, rnorm(totalN, -2, 1), rnorm(totalN, 2, 1))

  pop <- data.frame(id, X1, X2, X3, W1=W1[id], W2=W2[id])

  # ::::: 4) generate selection probabilities and potential outcome ::::::::

  E <- rnorm(sum(N), 0, 1)   # error terms for pot.

  pop$lps <- -0.5 + sel.coef*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$W2) # ps logit
  pop$Y0 <- 70 + trt.val*0 + out.coef*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$W2) + E
  pop$Y1 <- 70 + trt.val*1 + out.coef*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$W2) + E

  if (intX2W1) {
    pop$lps  <- pop$lps + sel.coef*pop$X2*pop$W1
    pop$Y0  <- pop$Y0 + out.coef*pop$X2*pop$W1
    pop$Y1  <- pop$Y1 + out.coef*pop$X2*pop$W1
  }

  if (intX1sq) {
    pop$lps  <- pop$lps + sel.coef*pop$X1^2
    pop$Y0  <- pop$Y0 + out.coef*pop$X1^2
    pop$Y1  <- pop$Y1 + out.coef*pop$X1^2
  }

  if (intX1X2) {
    pop$lps  <- pop$lps + sel.coef*pop$X1*pop$X2
    pop$Y0  <- pop$Y0 + out.coef*pop$X1*pop$X2
    pop$Y1  <- pop$Y1 + out.coef*pop$X1*pop$X2
  }

  if (intZW1) {
    pop$Y1  <- pop$Y1 + int.trt.val*pop$W1
  }

  if (intZX2) {
    pop$Y1  <- pop$Y1 + int.trt.val*pop$X2
  }

  if (intZX2sq) {
    pop$Y1  <- pop$Y1 + 0.5*pop$X2^2
  }

  if (intZX3) {
    pop$Y1  <- pop$Y1 + int.trt.val*pop$X3
  }

  pop$ps <- 1 / (1 + exp(-pop$lps))   # ps
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # treatment indicator
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))

  pop
}


ccrem.pop <- function(mu.w = c(0, 0), mu.j = c(0, 0),  mu.q = c(0, 0),
                      ids=1, Smpl.size = "40.50.30", trt.val=2, int.trt.val=1, sel.coef=0.3, out.coef=2,
                      intX2W1 = FALSE, intX1Q1 = FALSE, intX1M1 = FALSE, intX1X2 = FALSE, intZW1 = FALSE, intZX3 = FALSE, intX1sq=FALSE) {

  # ::::: 1) generate the number of observations in each level-2 Factor 1 :::::

  size <- strsplit(Smpl.size, "\\.")[[1]]  # level-2 Factor 1, level-2 Factor 2, level-1 units

  J <- as.numeric(size[1]) # level-2 Factor 1
  K <- as.numeric(size[2]) # level-2 Factor 2
  n.clus <- as.numeric(size[3]) # level-1 unit

  if (n.clus < 20) {  # variation of cluster sizes
    sd.clus <- 0.5
  } else {
    sd.clus <- 1
  }

  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes
  N[N<5] <- 5
  f1id <- factor(rep(ids:(J + ids - 1), N))             # level-2 Factor 1 id
  caseid <- factor(1:sum(N))

  pop <- data.frame(cbind(caseid, f1id))

  # ::::: 2) generate level-2 covariate, Ws (Factor 1) & Qs (Factor 2) :::::

  # Level-2 Factor 1 covariate, Ws
  W1 <- rnorm(J, 0, 1)
  W2 <- rnorm(J, 0, 1) # level-2 factor 1 unmeasured confounder
  names(W1) <- names(W2) <- as.factor(1:J)

  W.matrix <- cbind(1, W1, W2)

  # Level-2 Factor 2 covariate, Qs
  Q1 <- rnorm(K, 0, 1)
  Q2 <- rnorm(K, 0, 1) # level-2 factor 2 unmeasured confounder
  names(Q1) <- names(Q2) <- as.factor(1:K)

  Q.matrix <- cbind(1, Q1, Q2)

  # ::::: 3) create level-2 Factor 2 ID and cell-specfic ID :::::

  f1_order <- names(W2)
  f2_order <- names(Q2)

  for (i in 1:J) {

    if (i == 1) {

      temp.sm <- sample(temp.f1caseid <- pop[f1id==f1_order[i],1], round(0.6*N[as.numeric(f1_order[i])]))
      pop[temp.sm, "f2id"] <- f2_order[i]

      temp.sm3 <- temp.f1caseid[!(temp.f1caseid %in% temp.sm)]
      pop[temp.sm3, "f2id"] <- f2_order[i+1]


    } else if (i == J) {

      temp.sm <- sample(temp.f1caseid <- pop[f1id==f1_order[i],1], round(0.6*N[as.numeric(f1_order[i])]))
      pop[temp.sm, "f2id"] <- f2_order[i]

      temp.sm3 <- temp.f1caseid[!(temp.f1caseid %in% temp.sm)]
      pop[temp.sm3, "f2id"] <- f2_order[i-1]

    } else {

      temp.sm <- sample(temp.f1caseid <- pop[f1id==f1_order[i],1], round(0.5*N[as.numeric(f1_order[i])]))
      pop[temp.sm, "f2id"] <- f2_order[i]

      temp.sm2 <- sample(temp.f1caseid[!(temp.f1caseid %in% temp.sm)], round(0.25*N[as.numeric(f1_order[i])]))
      pop[temp.sm2, "f2id"] <- f2_order[i+1]

      temp.sm3 <- temp.f1caseid[!(temp.f1caseid %in% c(temp.sm, temp.sm2))]
      pop[temp.sm3, "f2id"] <- f2_order[i-1]

    }

  }

  pop$f2id <- factor(pop$f2id, levels= f2_order)
  pop$f1id <- factor(pop$f1id)
  pop$f12id = factor(with(pop, interaction(f1id, f2id, sep = "x")))
  f12id.num <- length(levels(pop$f12id))

  # ::::: 4) generate a level-1 covariate, X with Level2 Factor 1-specific means & covariates :::::

  totalN <- sum(N)

  X1 <- runif(totalN, -1, 1)
  X2 <- rnorm(totalN)
  X3_U <- runif(totalN, 0, 1)
  X3 <- ifelse(X3_U < 0.5, rnorm(totalN, -1, 1), rnorm(totalN, 1, 1))

  pop <- data.frame(pop, X1, X2, X3, W1=W1[f1id], W2=W2[f1id], Q1=Q1[pop$f2id], Q2=Q2[pop$f2id])

  # ::::: 5) generate selection probabilities and level-1 potential outcomes :::::

  E <- rnorm(totalN, 0, 1)   # error terms for pot.

  pop$lps <- -0.5 + sel.coef*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$W2 + pop$Q1 + pop$Q2) # ps logit
  pop$Y0 <- 70 + trt.val*0 + out.coef*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$W2 + pop$Q1 + pop$Q2) + E
  pop$Y1 <- 70 + trt.val*1 + out.coef*(pop$X1 + pop$X2 + pop$X3 + pop$W1 + pop$W2 + pop$Q1 + pop$Q2)  + E


  if (intX2W1) {
    pop$lps <- pop$lps + sel.coef*pop$X2*pop$W1
    pop$Y0 <- pop$Y0 + out.coef*pop$X2*pop$W1
    pop$Y1 <- pop$Y1 + out.coef*pop$X2*pop$W1
  }

  if (intX1Q1) {
    pop$lps <- pop$lps + sel.coef*pop$X1*pop$Q1
    pop$Y0 <- pop$Y0 + out.coef*pop$X1*pop$Q1
    pop$Y1 <- pop$Y1 + out.coef*pop$X1*pop$Q1
  }

  if (intX1sq) {
    pop$lps  <- pop$lps + sel.coef*pop$X1^2
    pop$Y0  <- pop$Y0 + out.coef*pop$X1^2
    pop$Y1  <- pop$Y1 + out.coef*pop$X1^2
  }

  if (intX1X2) {
    pop$lps <- pop$lps + sel.coef*pop$X1*pop$X2
    pop$Y0 <- pop$Y0 + out.coef*pop$X1*pop$X2
    pop$Y1 <- pop$Y1 + out.coef*pop$X1*pop$X2
  }

  if (intZW1) {
    pop$Y1 <- pop$Y1 + int.trt.val*pop$W1
  }

  if (intZX3) {
    pop$Y1 <- pop$Y1 + int.trt.val*pop$X3
  }


  # ::::: 6) generate actual selection and generate observed outcome :::::

  pop$ps <- 1 / (1 + exp(-pop$lps))   # ps
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # treatment indicator
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))

  pop$caseid <- NULL

  pop
}

twolevel_data = twolevel.pop(Smpl.size="20.50", intZW1 = T)
twolevel_data <- twolevel_data[, c("id", "Y", "Z", "X1",  "X2",  "X3",  "W1",  "W2", "lps", "ps", "Y0", "Y1")]
head(twolevel_data)

crossclassified_data = ccrem.pop(Smpl.size="20.20.100", intZW1 = T)
crossclassified_data <- crossclassified_data[, c("f1id",  "f2id",  "f12id", "Y", "Z", "X1",  "X2",  "X3",  "W1",  "W2", "Q1", "Q2", "lps", "ps", "Y0", "Y1")]
head(crossclassified_data)

usethis::use_data(twolevel_data, compress ="xz", overwrite = T)
usethis::use_data(crossclassified_data, compress ="xz", overwrite = T)
