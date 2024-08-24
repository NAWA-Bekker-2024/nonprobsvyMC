## R translations from https://github.com/cran/misclassGLM/blob/master/src/misclassGLMhelper.c

### loglikelihood
cflogit_validation <- function(par, data) {

  # Extract data elements
  y <- data$y
  x <- data$x
  setm <- data$setm
  w <- data$w
  px <- data$px
  ms <- data$ms
  pm <- data$pm

  # Calculate linear predictor
  beta <- par
  eta <- beta[1] + x %*% beta[2:(px+1)]

  # Add matrix component
  eta_mat <- eta + setm %*% beta[(px+2):(px+pm+1)]

  # Calculate probabilities
  p <- 1 / (1 + exp(-eta_mat))

  # Weighted probabilities
  wp <- w * p
  w1p <- w * (1-p)

  # Log likelihood components
  ll_y1 <- log(rowSums(wp))
  ll_y0 <- log(rowSums(w1p))

  # Full log likelihood
  ll <- sum(y * ll_y1 + (1-y) * ll_y0)

  # Return negative log likelihood
  return(-ll)
}

### gradient
cglogit_validation <- function(par, data) {

  # Extract data elements
  y <- data$y
  x <- data$x
  setm <- data$setm
  w <- data$w
  px <- data$px
  ms <- data$ms
  pm <- data$pm
  n <- length(y)
  p <- length(par)

  # Calculate linear predictor
  eta <- par[1] + x %*% par[2:(px+1)]

  # Add matrix component
  eta_mat <- sweep(setm %*% par[(px+2):(px+pm+1)], 1, eta, "+")

  # Calculate probabilities
  exp_eta <- exp(eta_mat)
  p_mat <- exp_eta / (1 + exp_eta)

  # Weighted probabilities
  wp <- w * p_mat
  w1p <- w * (1 - p_mat)

  # Log likelihood components
  lik <- rowSums(ifelse(y == 1, wp, w1p))

  # Calculate gradients
  tmp <- ifelse(y == 1, 1, -1) * exp_eta / (1 + exp_eta)^2 * w

  grad <- matrix(0, n, p)
  grad[, 1] <- rowSums(tmp)
  grad[, 2:(px+1)] <- x * rowSums(tmp)

  for (d in 1:pm) {
    grad[, px+1+d] <- rowSums(tmp * matrix(setm[d,], nrow=n, ncol=ms, byrow=TRUE))
  }

  # Adjust gradients by likelihood
  grad <- grad / lik

  # Sum gradients across observations
  ret <- -colSums(grad)

  return(ret)
}
