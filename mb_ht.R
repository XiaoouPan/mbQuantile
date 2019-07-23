## Goodness of fit testing via multiplier bootstrap
rm(list = ls())

library(quantreg)
library(MASS)
library(matrixStats)
library(logspline)
library(xtable)

multiBoot = function(X, Yadj, Yadj0, res, res0, n, alphaSeq, tau = 0.5, B = 1000) {
  mbStat = rep(0, B)
  for (b in 1:B) {
    w = rbinom(n, 1, 0.5)
    idx = which(w > 0)
    fit = rq(Yadj[idx] ~ X[idx, ], tau = tau)
    boot.res = as.numeric(fit$residuals)
    boot.res0 = Yadj0[idx] - median(Yadj0[idx])
    mbStat[b] = 2 * (sum(abs(boot.res0)) - sum(abs(boot.res))) - 
      2 * (sum(abs(res0[idx])) - sum(abs(res[idx])))
  }
  return (quantile(mbStat, probs = 1 - alphaSeq))
}

geneMulti = function(n, d) {
  V = runif(d, 0.5, 1)
  Sigma = diag(V)
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      Sigma[i, j] = Sigma[j, i] = 0.5^(j - i) * sqrt(V[i] * V[j])
    }
  }
  return (mvrnorm(n, rep(0, d), Sigma))
}

geneEquCor = function(n, d) {
  V = runif(d, 0.5, 1)
  Sigma = diag(d)
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      Sigma[i, j] = Sigma[j, i] = 0.5 * sqrt(V[i] * V[j])
    }
  }
  return (mvrnorm(n, rep(0, d), Sigma))
}

geneMixNoise = function(n) {
  e1 = rnorm(n, -1, 1)
  e2 = rnorm(n, 1, 1)
  a = rbinom(n, 1, 0.5)
  return (a * e1 + (1 - a) * e2)
}

getContNoise = function(n) {
  e1 = rnorm(n)
  e2 = rnorm(n, 0, 10)
  a = rbinom(n, 1, 0.95)
  return (a * e1 + (1 - a) * e2)
}

n = 200
d = 15
B = 1000
M = 200
tau = 0.5
alphaSeq = c(0.01, 0.05, 0.1)

#model = "homo"
model = "hetero"

rej = matrix(0, 3, 12)
pb = txtProgressBar(style = 3)
for (m in 1:M) {
  set.seed(m)
  # Independent X
  X = matrix(rnorm(n * d), n, d)
  # Multivariate X with decreasing correlation
  #X = geneMulti(n, d)
  # Multivariate X with equal correlation
  #X = geneEquCor(n, d)
  
  # t error
  error = rt(n, df = 2)
  # mixture normal error
  #error = geneMixNoise(n)
  # contaminated error
  #error = getContNoise(n)
  
  ## Null setting
  betaStar = rep(0, d + 1)
  Y = NULL
  if (model == "homo") {
    Y = as.numeric(cbind(rep(1, n), X) %*% betaStar) + error
  } else if (model == "hetero") {
    gamma = 1
    D = 2 * diag(exp(gamma * X[, 1]) / (1 + exp(gamma * X[, 1])))
    Y = as.numeric(cbind(rep(1, n), X) %*% betaStar) + as.numeric(D %*% error) 
  }
  fit = rq(Y ~ X, tau = tau)
  res = as.numeric(fit$residuals)
  #adj = hat(X) * (tau - (res < 0)) / akj(res, z = 0)$dens
  #Yadj = Y - adj
  res0 = Y - median(Y)
  stat = sum(abs(res0)) - sum(abs(res))
  fit0 = rq(Y ~ 1, tau = tau)
  #adj0 = (1 / n) * (tau - (res0 < 0)) / akj(res0, z = 0)$dens
  #Yadj0 = Y - adj0
  ## 1. ANOVA, Wald test
  test = anova(fit, fit0, test = "Wald")
  rej[, 1] = rej[, 1] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 2. ANOVA, Gutenbrunner etal
  test = anova(fit, fit0, test = "rank")
  rej[, 2] = rej[, 2] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 3. ANOVA, Chen etal
  test = anova(fit, fit0, test = "anowar", R = B)
  rej[, 3] = rej[, 3] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 4. Multiplier bootstrap
  criVal = multiBoot(X, Y, Y, res, res0, n, alphaSeq, tau, B)
  rej[, 4] = rej[, 4] + as.numeric(stat > as.numeric(criVal))
  
  ## Sparse setting
  s = 1
  betaStar = c(0, rep(0.5, s), rep(0, d - s))
  Y = NULL
  if (model == "homo") {
    Y = as.numeric(cbind(rep(1, n), X) %*% betaStar) + error
  } else if (model == "hetero") {
    gamma = 1
    D = 2 * diag(exp(gamma * X[, 1]) / (1 + exp(gamma * X[, 1])))
    Y = as.numeric(cbind(rep(1, n), X) %*% betaStar) + as.numeric(D %*% error)
  }
  fit = rq(Y ~ X, tau = tau)
  res = as.numeric(fit$residuals)
  #f0 = akj(res, z = 0)$dens
  #adj = hat(X) * (tau - (res < 0)) / f0
  #Yadj = Y - adj
  res0 = Y - median(Y)
  stat = sum(abs(res0)) - sum(abs(res))
  fit0 = rq(Y ~ 1, tau = tau)
  #adj0 = (1 / n) * (tau - (res0 < 0)) / akj(res0, z = 0)$dens
  #Yadj0 = Y - adj0
  ## 1. ANOVA, Wald test
  test = anova(fit, fit0, test = "Wald")
  rej[, 5] = rej[, 5] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 2. ANOVA, Gutenbrunner etal
  test = anova(fit, fit0, test = "rank")
  rej[, 6] = rej[, 6] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 3. ANOVA, Chen etal
  test = anova(fit, fit0, test = "anowar", R = B)
  rej[, 7] = rej[, 7] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 4. Multiplier bootstrap
  criVal = multiBoot(X, Y, Y, res, res0, n, alphaSeq, tau, B)
  rej[, 8] = rej[, 8] + as.numeric(stat > as.numeric(criVal))
  
  ## dense setting
  s = 10
  betaStar = c(0, rep(0.1, s), rep(0, d - s))
  Y = NULL
  if (model == "homo") {
    Y = as.numeric(cbind(rep(1, n), X) %*% betaStar) + error
  } else if (model == "hetero") {
    gamma = 1
    D = 2 * diag(exp(gamma * X[, 1]) / (1 + exp(gamma * X[, 1])))
    Y = as.numeric(cbind(rep(1, n), X) %*% betaStar) + as.numeric(D %*% error)
  }
  fit = rq(Y ~ X, tau = tau)
  res = as.numeric(fit$residuals)
  #f0 = akj(res, z = 0)$dens
  #adj = hat(X) * (tau - (res < 0)) / f0
  #Yadj = Y - adj
  res0 = Y - median(Y)
  stat = sum(abs(res0)) - sum(abs(res))
  fit0 = rq(Y ~ 1, tau = tau)
  #adj0 = (1 / n) * (tau - (res0 < 0)) / akj(res0, z = 0)$dens
  #Yadj0 = Y - adj0
  ## 1. ANOVA, Wald test
  test = anova(fit, fit0, test = "Wald")
  rej[, 9] = rej[, 9] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 2. ANOVA, Gutenbrunner etal
  test = anova(fit, fit0, test = "rank")
  rej[, 10] = rej[, 10] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 3. ANOVA, Chen etal
  test = anova(fit, fit0, test = "anowar", R = B)
  rej[, 11] = rej[, 11] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
  ## 4. Multiplier bootstrap
  criVal = multiBoot(X, Y, Y, res, res0, n, alphaSeq, tau, B)
  rej[, 12] = rej[, 12] + as.numeric(stat > as.numeric(criVal))
  setTxtProgressBar(pb, m / M)
}

rej = data.frame(rej / M)
names(rej) = rep(c("Wald", "rank", "mb-exp", "mb-Rad"), 3)
row.names(rej) = c("0.01", "0.05", "0.1")
rej
