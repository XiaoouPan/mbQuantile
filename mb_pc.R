##########################################
#### Simulation codes for power curve ####
##########################################
rm(list = ls())

library(quantreg)
library(MASS)
library(matrixStats)
library(logspline)
library(xtable)
library(tikzDevice)

multiBoot = function(X, Yadj, res, res0, n, alphaSeq, tau = 0.5, B = 1000) {
  mbStat = rep(0, B)
  for (b in 1:B) {
    w = rbinom(n, 1, 0.5)
    idx = which(w > 0)
    fit = rq(Yadj[idx] ~ X[idx, ], tau = tau)
    boot.res = as.numeric(fit$residuals)
    boot.res0 = Yadj[idx] - median(Yadj[idx])
    mbStat[b] = 2 * (sum(abs(boot.res0)) - sum(abs(boot.res))) - 2 * (sum(abs(res0[idx])) - sum(abs(res[idx])))
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
  Sigma = diag(V)
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
  e2 = rnorm(n, 0, 5)
  a = rbinom(n, 1, 0.9)
  return (a * e1 + (1 - a) * e2)
}

n = 200
d = 15
B = 1000
M = 200
tau = 0.5
alphaSeq = 0.05
len = 20

model = "homo"
#model = "hetero"

## Homo, mn1, indep, sparse setting
s = 1
signal = seq(0.1, 1, length.out = len)
## Homo, mn1, indep, dense setting
#s = 10
#signal = seq(0.001, 0.35, length.out = len)
## Hetero, mn1, indep, sparse setting
#s = 1
#signal = seq(0.001, 0.8, length.out = len)
## Hetero, mn1, indep, dense setting
#s = 10
#signal = seq(0.001, 0.35, length.out = len)

pw1 = matrix(0, M, len)
pw2 = matrix(0, M, len)
pw3 = matrix(0, M, len)
pb = txtProgressBar(style = 3)
for (m in 1:M) {
  set.seed(m)
  # Independent X
  X = matrix(rnorm(n * d), n, d)
  # type I mixture normal error
  error = geneMixNoise(n)

  for (i in 1:len) {
    betaStar = c(0, rep(signal[i], s), rep(0, d - s))
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
    res0 = Y - median(Y)
    stat = sum(abs(res0)) - sum(abs(res))
    fit0 = rq(Y ~ 1, tau = tau)
    ## 1. Gutenbrunner etal
    test = anova(fit, fit0, test = "rank")
    pw1[m, i] = pw1[m, i] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
    ## 2. Chen etal
    test = anova(fit, fit0, test = "anowar", R = B)
    pw2[m, i] = pw2[m, i] + as.numeric(as.numeric(test$table[4]) < alphaSeq)
    ## 3. Multiplier bootstrap
    criVal = multiBoot(X, Y, res, res0, n, alphaSeq, tau, B)
    pw3[m, i] = pw3[m, i] + as.numeric(stat > as.numeric(criVal))
    
    setTxtProgressBar(pb, ((m - 1) * len + i) / (M * len))
  }
}

rst = cbind(pw1, pw2, pw3)
write.csv(rst, "Results/hetero_mn1_indep_den_n350.csv")

pc1 = colMeans(pw1)
pc2 = colMeans(pw2)
pc3 = colMeans(pw3)
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
plot(signal, pc1, type = "l", lwd = 3, xlim = range(signal), ylim = c(0, 1), xlab = "Signal strength", ylab = "Power")
lines(signal, pc2, col = "blue", lwd = 3)
lines(signal, pc3, col = "red", lwd = 3)
colors = c("black", "blue", "red")
labels = c("rank", "mb-exp", "mb-Rad")
legend("bottomright", labels, col = colors, lwd = 3, cex = 1.8)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)
