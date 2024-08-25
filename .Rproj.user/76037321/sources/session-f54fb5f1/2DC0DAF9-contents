cfgauss_validation <- function(beta, data) {
  y <- data$y
  x <- data$x
  setm <- data$setm
  w <- data$w
  n <- data$n
  px <- data$px
  ms <- data$ms
  pm <- data$pm

  sigma2 <- -2 * beta[px+pm+2]^2
  lnsqrt2pi <- log(sqrt(2*pi) * beta[px+pm+2])

  eta <- beta[1] + x %*% beta[2:(px+1)]
  eta_mat <- sweep(setm %*% beta[(px+2):(px+pm+1)], 1, eta, "+")

  tmp <- rowSums(exp((y - eta_mat)^2 / sigma2) * w)

  ret <- sum(log(tmp)) - n*lnsqrt2pi
  return(-ret)
}

cggauss_validation <- function(beta, data) {
  y <- data$y
  x <- data$x
  setm <- data$setm
  w <- data$w
  n <- data$n
  px <- data$px
  ms <- data$ms
  pm <- data$pm

  sigma2a <- beta[px+pm+2]^2
  sigma2 <- -2 * sigma2a

  eta <- beta[1] + x %*% beta[2:(px+1)]
  eta_mat <- sweep(setm %*% beta[(px+2):(px+pm+1)], 1, eta, "+")

  predicts <- y - eta_mat
  exppredicts <- exp(predicts^2 / sigma2) * w

  tempsumm <- rowSums(exppredicts)
  tmp <- rowSums(predicts * exppredicts) / tempsumm

  ret <- c(sum(tmp),
           colSums(tmp * x),
           colSums((predicts * exppredicts) %*% t(setm)) / tempsumm,
           sum((rowSums(predicts^2 * exppredicts) / sigma2a - 1) / tempsumm))

  ret[1:(px+pm+1)] <- ret[1:(px+pm+1)] / sigma2a
  ret[px+pm+2] <- ret[px+pm+2] / beta[px+pm+2]

  return(-ret)
}

cflogit_validation <- function(beta, data) {
  y <- data$y
  x <- data$x
  setm <- data$setm
  w <- data$w
  n <- data$n
  px <- data$px
  ms <- data$ms
  pm <- data$pm

  eta <- beta[1] + x %*% beta[2:(px+1)]
  eta_mat <- sweep(setm %*% beta[(px+2):(px+pm+1)], 1, eta, "+")

  p <- 1 / (1 + exp(-eta_mat))

  tmp <- rowSums(ifelse(y == 1, p * w, (1-p) * w))
  ret <- sum(log(tmp))

  return(-ret)
}

cglogit_validation <- function(beta, data) {
  y <- data$y
  x <- data$x
  setm <- data$setm
  w <- data$w
  n <- data$n
  px <- data$px
  ms <- data$ms
  pm <- data$pm

  eta <- beta[1] + x %*% beta[2:(px+1)]
  eta_mat <- sweep(setm %*% beta[(px+2):(px+pm+1)], 1, eta, "+")

  p <- 1 / (1 + exp(-eta_mat))

  lik <- rowSums(ifelse(y == 1, p * w, (1-p) * w))

  tmp <- ifelse(y == 1, 1, -1) * p * (1-p) * w

  ret <- c(sum(rowSums(tmp) / lik),
           colSums((rowSums(tmp) / lik) * x),
           colSums((tmp / lik) %*% t(setm)))

  return(-ret)
}
