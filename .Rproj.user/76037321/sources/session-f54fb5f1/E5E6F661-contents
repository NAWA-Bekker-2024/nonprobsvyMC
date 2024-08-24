#' Function to correct misclassification based on `nonprobsvy` object
#' @import
#'
#' @param obj `nonprobsvy` object (class: `nonprobsvy`)
#' @param y either a formula with variables that are misclassified or a list of objects if validation sample is provided
#' @param x either a formula with variables that are misclassified or a list of objects if validation sample is provided
#' @param y_class_mat either a single matrix or a list of misclassification matrices for `y` variables
#' @param x_class_mat either a single matrix or a list of misclassification matrices for `x` variables
#' @param valid_sample a data.frame with validation sample where columns are determined by `y` and `x` arguments
#' @param method a vector of possible methods to correct for the misclassification
#' @param svy a `survey.design2` object
#'
#' @return An `nonprobsvy` object corrected for misclassification
#'
#' @details
#' Details will be provided here
#'
#' @examples
#' # examples will be here
#'
#' @export
correct_misclass <- function(obj,
                             y,
                             x,
                             y_class_mat,
                             x_class_mat,
                             valid_sample,
                             svy, ## currently need to be provided
                             method = "simex",
                             control = list(simex.lambda = c(0.5, 1, 1.5, 2),
                                            simex.boot = 50)) {

  split_formula <- function(formula) {
    formula_char <- deparse(formula)
    parts <- strsplit(formula_char, "~")[[1]]
    rhs <- trimws(parts[2])
    lhs <- strsplit(trimws(parts[1]), "\\+")[[1]]
    lhs <- trimws(lhs)
    formula_list <- list()
    for (response in lhs) {
      new_formula <- as.formula(paste(response, "~", rhs))
      formula_list[[response]] <- new_formula
    }
    return(formula_list)
  }

  obj_class <- class(obj)
  if (obj_class[1] != "nonprobsvy") {
    stop("Input must be object of `nonprobsvy` class")
  }
  obj_method <- gsub("nonprobsvy_", "", obj_class[2])
  stopifnot("Currently, only mass imputation is supported" = obj_method == "mi")

  ## processing variables
  ys <- obj$y  # List of Y variables
  glm_objs <- obj$outcome  # List of GLM objects
  glm_objs_formulas <- split_formula(obj$call$outcome)

  results <- list()

  for (y_name in names(ys)) {
    glm_obj <- glm_objs[[y_name]]

    # Update GLM object
    glm_obj$call$formula <- glm_objs_formulas[[y_name]]
    glm_obj$call$family <- obj$call$family_outcome
    glm_obj$call$control <- list(obj$control$control_outcome$epsilon,
                                 obj$control$control_outcome$maxit,
                                 obj$control$control_outcome$trace)

    glm_obj$call$start <- coef(glm_obj)
    glm_obj_data <- glm_obj$data

    # SIMEX procedure
    simex_lin <- simex::mcsimex(model = glm_obj,
                                SIMEXvariable = names(x_class_mat),
                                mc.matrix = x_class_mat,
                                asymptotic = FALSE,
                                jackknife.estimation = FALSE,
                                fitting.method = "linear",
                                B = control$simex.boot,
                                lambda = control$simex.lambda)

    simex_log <- simex::refit(simex_lin, fitting.method = "loglinear",
                              asymptotic = FALSE, jackknife.estimation = FALSE)
    simex_qua <- simex::refit(simex_lin, fitting.method = "quadratic",
                              asymptotic = FALSE, jackknife.estimation = FALSE)

    # Update survey object
    svy$variables[[paste0(y_name, "_lin")]] <- predict(simex_lin, newdata=svy$variables, type = "response")
    svy$variables[[paste0(y_name, "_log")]] <- predict(simex_log, newdata=svy$variables, type = "response")
    svy$variables[[paste0(y_name, "_qua")]] <- predict(simex_qua, newdata=svy$variables, type = "response")

    # Calculate means and confidence intervals
    formula_str <- paste0("~", y_name, "_lin+", y_name, "_log+", y_name, "_qua")
    mu_mi <- svymean(as.formula(formula_str), svy)
    mu_mi_ci <- confint(mu_mi)

    result <- data.frame(mu_mi)
    result$ci_low <- mu_mi_ci[,1]
    result$ci_upp <- mu_mi_ci[,2]

    results[[y_name]] <- result
  }

  names(results) <- NULL
  results <- do.call("rbind", results)
  results$Y <- rownames(results)
  rownames(results) <- NULL

  return(results[, c("Y", "mean", "SE", "ci_low", "ci_upp")])
}
