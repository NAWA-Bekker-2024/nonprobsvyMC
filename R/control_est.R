#' Control parameters for estimation
#'
#' This function creates a list of control parameters for use in estimation
#' procedures. It provides defaults based on the `stats::optim` function and
#' includes additional parameters from `simex::mcsimex`.
#'
#' @param optim.maxit Maximum number of iterations for optim. Default is 100.
#' @param optim.abstol Absolute convergence tolerance for optim. Default is sqrt(.Machine$double.eps).
#' @param optim.reltol Relative convergence tolerance for optim. Default is sqrt(.Machine$double.eps).
#' @param optim.trace Non-negative integer. If positive, tracing information is produced.
#'                    Higher values may produce more tracing information. Default is 0.
#' @param simex.lambda A vector of λ values used in the simulation step. Default is c(0.5, 1, 1.5, 2).
#' @param simex.B Number of simulation runs. Default is 100.
#' @param simex.fitting.method Method used for fitting. Default is "quadratic".
#' @param simex.asymptotic Logical indicating whether asymptotic variance should be calculated. Default is FALSE.
#' @param simex.variance.jackknife Logical indicating whether jackknife variance estimation should be used. Default is FALSE.
#'
#' @return A list with the specified control parameters.
#'
#' @examples
#' ctrl <- control_est(optim.maxit = 200, optim.trace = 1, simex.lambda = c(0.5, 1, 1.5, 2, 2.5), simex.B = 200)
#'
#' @export
control_est <- function(optim.maxit = 100,
                        optim.abstol = sqrt(.Machine$double.eps),
                        optim.reltol = sqrt(.Machine$double.eps),
                        optim.trace = 0,
                        optim.fnscale = 1,
                        simex.lambda = c(0.5, 1, 1.5, 2),
                        simex.B = 100,
                        simex.fitting.method = "quadratic",
                        simex.asymptotic = FALSE,
                        simex.variance.jackknife = FALSE) {
  list(optim.maxit = optim.maxit,
       optim.abstol = optim.abstol,
       optim.reltol = optim.reltol,
       optim.trace = optim.trace,
       simex.lambda = simex.lambda,
       simex.B = simex.B,
       simex.fitting.method = simex.fitting.method,
       simex.asymptotic = simex.asymptotic,
       simex.variance.jackknife = simex.variance.jackknife)
}
