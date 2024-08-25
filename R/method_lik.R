#' Internal functions correct misclassification using likelihood methods (one variable only)
#'
#' @param y target variable
#' @param x auxiliary variables (measured without error)
#' @param mat_z auxiliary variable measured with error (single)
#' @param mat_p matrix with probabilities for each level of z
#' @param weights frequency weights for each case
#' @param family target variable family
#' @param method non-probability method (either mass imputation, inverse probability weighting or doubly robust estimators)
#' @param start starting parameters
#' @param control control parameters passed to optim
#'
#' @example
#' ## example here
method_lik <- function(y, x, mat_z, mat_p, weights, family, method, start, control) {

  ## binomial distribution
  mi_binom_ll <- function(par, y, x, mat_z, mat_p, weights) {

    betas1 <- par[1:ncol(x)]
    betas2 <- par[(ncol(x)+1):length(par)]
    eta <- as.numeric(x %*% betas1)
    eta_mat <- rep(1, NROW(eta)) %*% (betas2 %*% t(mat_z)) + eta
    probs <- plogis(eta_mat)
    tmp <- rowSums((probs*y + (1-y)*(1-probs))*mat_p*weights)
    return(-sum(log(tmp)))
  }

  mi_binom_gr <- function(par, y, x, mat_z, mat_p, weights) {
    betas1 <- par[1:ncol(x)]
    betas2 <- par[(ncol(x)+1):length(par)]
    eta <- as.numeric(x %*% betas1)
    eta_mat <- rep(1, NROW(eta)) %*% (betas2 %*% t(mat_z)) + eta
    probs <- plogis(eta_mat)
    probs <- plogis(eta_mat)
    lik <- rowSums((probs*y + (1-y)*(1-probs))*mat_p*weights)

    tmp <- ifelse(y == 1, 1, -1) * probs * (1-probs) * mat_p

    ret <- c(sum(rowSums(tmp) / lik),
             colSums((rowSums(tmp) / lik) * x),
             colSums((tmp / lik) %*% setm))

    grad <- c(colSums((rowSums(tmp) / lik)*x*weights),
              colSums(((tmp / lik) %*% mat_z)*weights))

    return(-grad)
  }

  ## gaussian distribution (for verification)

  mi_gaussian_ll <- function(par, y, x, mat_z, mat_p, weights) {
    betas1 <- par[1:(ncol(x) + 1)]  # +1 for sigma
    betas2 <- par[(ncol(x) + 2):length(par)]
    sigma <- exp(betas1[length(betas1)])  # Ensure sigma is positive
    betas1 <- betas1[-length(betas1)]  # Remove sigma from betas1

    eta <- as.numeric(x %*% betas1)
    eta_mat <- rep(1, NROW(eta)) %*% (betas2 %*% t(mat_z)) + eta

    tmp <- rowSums(dnorm(y, mean = eta_mat, sd = sigma) * mat_p * weights)
    return(-sum(log(tmp)))
  }

  mi_gaussian_gr <- function(par, y, x, mat_z, mat_p, weights) {
    betas1 <- par[1:(ncol(x) + 1)]  # +1 for sigma
    betas2 <- par[(ncol(x) + 2):length(par)]
    sigma <- exp(betas1[length(betas1)])  # Ensure sigma is positive
    betas1 <- betas1[-length(betas1)]  # Remove sigma from betas1

    eta <- as.numeric(x %*% betas1)
    eta_mat <- rep(1, NROW(eta)) %*% (betas2 %*% t(mat_z)) + eta

    residuals <- y - eta_mat
    lik <- rowSums(dnorm(y, mean = eta_mat, sd = sigma) * mat_p * weights)

    tmp <- (residuals / sigma^2) * mat_p

    grad_betas1 <- colSums((rowSums(tmp) / lik) * x * weights)
    grad_betas2 <- colSums(((tmp / lik) %*% mat_z) * weights)

    grad_sigma <- sum(((rowSums((residuals^2 / sigma^3 - 1 / sigma) * mat_p) / lik) * weights))

    grad <- c(grad_betas1, grad_sigma * sigma, grad_betas2)  # Chain rule for sigma

    return(-grad)
  }

  ## poisson distribution (for verification)
  mi_poisson_ll <- function(par, y, x, mat_z, mat_p, weights) {
    betas1 <- par[1:ncol(x)]
    betas2 <- par[(ncol(x)+1):length(par)]
    eta <- as.numeric(x %*% betas1)
    eta_mat <- rep(1, NROW(eta)) %*% (betas2 %*% t(mat_z)) + eta
    lambda <- exp(eta_mat)
    tmp <- rowSums((dpois(y, lambda) * mat_p * weights))
    return(-sum(log(tmp)))
  }

  mi_poisson_gr <- function(par, y, x, mat_z, mat_p, weights) {
    betas1 <- par[1:ncol(x)]
    betas2 <- par[(ncol(x)+1):length(par)]
    eta <- as.numeric(x %*% betas1)
    eta_mat <- rep(1, NROW(eta)) %*% (betas2 %*% t(mat_z)) + eta
    lambda <- exp(eta_mat)

    lik <- rowSums(dpois(y, lambda) * mat_p * weights)
    tmp <- (y - lambda) * mat_p

    grad <- c(
      colSums((rowSums(tmp) / lik) * x * weights),
      colSums(((tmp / lik) %*% mat_z) * weights)
    )

    return(-grad)
  }

  ## here are switches for family but this will be also required for ipw and dr
  opt_fun <- switch(family()$family,
                    binomial = mi_binom_ll,
                    gaussian = mi_gaussian_ll,
                    poissomn = mi_poisson_ll)

  opt_gr <- switch(family()$family,
                   binomial = mi_binom_gr,
                   gaussian = mi_gaussian_gr,
                   poissomn = mi_poisson_gr)

  ## optim with gradient here
  optim_res <- stats::optim(par = start,
                            fn = opt_fun,
                            gr = opt_gr,
                            method = "BFGS",
                            hessian = T,
                            control = list(maxit = control$optim.maxit,
                                           abstol = control$optim.abstol,
                                           reltol = control$optim.reltol,
                                           trace = control$optim.trace),
                            ## parameters for opt_fun and opt_gr
                            y=y, x=x, mat_z=mat_z, mat_p=mat_p, weights=weights
                            )
  return(optim_res)
}
