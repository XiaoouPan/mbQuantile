## Confidence interval via multiplier bootstrap
rm(list = ls())

library(quantreg)
library(MASS)
library(matrixStats)
library(xtable)

getPercCI = function(estBoot, alphaSeq) {
  rst = NULL
  for (alpha in alphaSeq) {
    q1 = rowQuantiles(estBoot, probs = alpha / 2)
    q2 = rowQuantiles(estBoot, probs = 1 - alpha / 2)
    rst = cbind(rst, q1, q2)
  }
  return (rst) 
}

getPivCI = function(est, estBoot, alphaSeq) {
  rst = NULL
  for (alpha in alphaSeq) {
    q1 = rowQuantiles(estBoot, probs = alpha / 2)
    q2 = rowQuantiles(estBoot, probs = 1 - alpha / 2)
    rst = cbind(rst, 2 * est - q2, 2 * est - q1)
  }
  return (rst) 
}

getNormCI = function(est, sd, alphaSeq) {
  rst = NULL
  for (alpha in alphaSeq) {
    z = qnorm(1 - alpha / 2)
    lower = est - z * sd
    upper = est + z * sd
    rst = cbind(rst, lower, upper)
  }
  return (rst)
}

pairBoot = function(fit, betaHat, alphaSeq, B = 1000) {
  start = proc.time()
  list = summary(fit, se = "boot", bsmethod = "xy", R = B)
  se = as.numeric(list$coefficients[-1, 2])
  CI = getNormCI(betaHat, se, alphaSeq)
  end = proc.time()
  elapsed = as.numeric(end - start)[3]
  return (list(CI = CI, elapsed = elapsed))
}

efBoot = function(fit, betaHat, alphaSeq, B = 1000) {
  start = proc.time()
  list = summary(fit, se = "boot", bsmethod = "pwy", R = B)
  se = as.numeric(list$coefficients[-1, 2])
  CI = getNormCI(betaHat, se, alphaSeq)
  end = proc.time()
  elapsed = as.numeric(end - start)[3]
  return (list(CI = CI, elapsed = elapsed))
}

wildBoot = function(fit, betaHat, alphaSeq, B = 1000) {
  start = proc.time()
  list = summary(fit, se = "boot", bsmethod = "wild", R = B)
  se = as.numeric(list$coefficients[-1, 2])
  CI = getNormCI(betaHat, se, alphaSeq)
  end = proc.time()
  elapsed = as.numeric(end - start)[3]
  return (list(CI = CI, elapsed = elapsed))
}

multiBoot = function(X, Y, betaHat, n, d, alphaSeq, tau = 0.5, B = 1000) {
  multiBeta = matrix(0, d, B)
  start = proc.time()
  for (b in 1:B) {
    w = rbinom(n, 1, 0.5)
    idx = which(w > 0)
    fit = rq(Y[idx] ~ X[idx, ], tau = tau)
    multiBeta[, b] = as.numeric(fit$coefficients)[-1]
  }
  percCI = getPercCI(multiBeta, alphaSeq)
  normCI = getNormCI(betaHat, rowSds(multiBeta), alphaSeq)
  end = proc.time()
  elapsed = as.numeric(end - start)[3]
  return (list(percCI = percCI, normCI = normCI, elapsed = elapsed))
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
  e2 = rnorm(n, 1, 5)
  a = rbinom(n, 1, 0.9)
  return (a * e1 + (1 - a) * e2)
}

getSummary = function(lower, upper, elipsed, d, M) {
  rst = matrix(0, d * 5, 2)
  rst[, 1] = rowMeans(lower <= rep(betaStar[-1], 5) & upper >= rep(betaStar[-1], 5))
  rst[, 2] = rowMeans(upper - lower)
  rst = matrix(rst, d, 2 * 5)
  elapsed = elapsed / M
  return (list(rst = rst, elapsed = elapsed))
}

n = 200
d = 10
B = 1000
M = 200
tau = 0.5
betaStar = rep(2, d + 1)

alphaSeq = c(0.05, 0.1, 0.2)

lower005 = matrix(0, 5 * d, M)
upper005 = matrix(0, 5 * d, M)
lower01 = matrix(0, 5 * d, M)
upper01 = matrix(0, 5 * d, M)
lower02 = matrix(0, 5 * d, M)
upper02 = matrix(0, 5 * d, M)
elapsed = rep(0, 5)

pb = txtProgressBar(style = 3)
for (m in 1:M) {
  set.seed(m)
  # Independent X
  #X = matrix(rnorm(n * d), n, d)
  # Multivariate X with decreasing correlation
  #X = geneMulti(n, d)
  # Multivariate X with equal correlation
  X = geneEquCor(n, d)
  
  # t error
  #error = rt(n, df = 2)
  # mixture normal error
  #error = geneMixNoise(n)
  # contaminated error
  error = getContNoise(n)
  
  # homo model
  #Y = as.numeric(cbind(rep(1, n), X) %*% betaStar) + error
  # hetero model
  gamma = 1
  D = 2 * diag(exp(gamma * X[, 1]) / (1 + exp(gamma * X[, 1])))
  Y = as.numeric(cbind(rep(1, n), X) %*% betaStar) + as.numeric(D %*% error)
  
  fit = rq(Y ~ X, tau = tau)
  betaHat = as.numeric(fit$coefficients)[-1]
  res = as.numeric(fit$residuals)
  f0 = akj(res, z = 0)$dens
  adj = hat(X) * (tau - (res < 0)) / f0
  Yadj = Y - adj
  ## 1. Pairwise bootstrap, ET
  list = pairBoot(fit, betaHat, alphaSeq, B)
  lower005[1:d, m] = list$CI[, 1]
  upper005[1:d, m] = list$CI[, 2]
  lower01[1:d, m] = list$CI[, 3]
  upper01[1:d, m] = list$CI[, 4]
  lower02[1:d, m] = list$CI[, 5]
  upper02[1:d, m] = list$CI[, 6]
  elapsed[1] = elapsed[1] + list$elapsed 
  ## 2. Estimation function bootstrap, PWY
  list = efBoot(fit, betaHat, alphaSeq, B)
  lower005[(d + 1):(2 * d), m] = list$CI[, 1]
  upper005[(d + 1):(2 * d), m] = list$CI[, 2]
  lower01[(d + 1):(2 * d), m] = list$CI[, 3]
  upper01[(d + 1):(2 * d), m] = list$CI[, 4]
  lower02[(d + 1):(2 * d), m] = list$CI[, 5]
  upper02[(d + 1):(2 * d), m] = list$CI[, 6]
  elapsed[2] = elapsed[2] + list$elapsed 
  ## 3. Wild bootstrap, FHH
  list = wildBoot(fit, betaHat, alphaSeq, B)
  lower005[(2 * d + 1):(3 * d), m] = list$CI[, 1]
  upper005[(2 * d + 1):(3 * d), m] = list$CI[, 2]
  lower01[(2 * d + 1):(3 * d), m] = list$CI[, 3]
  upper01[(2 * d + 1):(3 * d), m] = list$CI[, 4]
  lower02[(2 * d + 1):(3 * d), m] = list$CI[, 5]
  upper02[(2 * d + 1):(3 * d), m] = list$CI[, 6]
  elapsed[3] = elapsed[3] + list$elapsed
  ## 4. Multiplier bootstrap, percentile CI
  list = multiBoot(X, Yadj, betaHat, n, d, alphaSeq, tau, B)
  lower005[(3 * d + 1):(4 * d), m] = list$percCI[, 1]
  upper005[(3 * d + 1):(4 * d), m] = list$percCI[, 2]
  lower01[(3 * d + 1):(4 * d), m] = list$percCI[, 3]
  upper01[(3 * d + 1):(4 * d), m] = list$percCI[, 4]
  lower02[(3 * d + 1):(4 * d), m] = list$percCI[, 5]
  upper02[(3 * d + 1):(4 * d), m] = list$percCI[, 6]
  elapsed[4] = elapsed[4] + list$elapsed
  ## 5. Multiplier bootstrap, normal CI
  lower005[(4 * d + 1):(5 * d), m] = list$normCI[, 1]
  upper005[(4 * d + 1):(5 * d), m] = list$normCI[, 2]
  lower01[(4 * d + 1):(5 * d), m] = list$normCI[, 3]
  upper01[(4 * d + 1):(5 * d), m] = list$normCI[, 4]
  lower02[(4 * d + 1):(5 * d), m] = list$normCI[, 5]
  upper02[(4 * d + 1):(5 * d), m] = list$normCI[, 6]
  elapsed[5] = elapsed[5] + list$elapsed
  setTxtProgressBar(pb, m / M)
}

report = NULL
sum005 = getSummary(lower005, upper005, elapsed, d, M)$rst
report = rbind(report, colMeans(sum005))
sum01 = getSummary(lower01, upper01, elapsed, d, M)$rst
report = rbind(report, colMeans(sum01))
sum02 = getSummary(lower02, upper02, elapsed, d, M)$rst
report = rbind(report, colMeans(sum02))
report = as.data.frame(report)
names(report) = rep(c("pair", "pwy", "wild", "mb-per", "mb-norm"), 2)
row.names(report) = c("0.05: ", "0.1: ", "0.2: ")
report

