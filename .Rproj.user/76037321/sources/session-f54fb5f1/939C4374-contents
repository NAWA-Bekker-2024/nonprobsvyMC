## function for optimization
##
library(nonprobsvy)
library(nonprobsvyMC)
library(rootSolve)

## without classsification errors
mle_logistic <- function(beta, x_a, x_b, d, p_mat) {
  pr_a <- as.numeric(exp(x_a %*% beta))
  pr_b <- as.numeric(exp(x_b %*% beta))
  ll <- sum(log(pr_a)) - sum(d*log(1 + pr_b))
  return(-ll)
}

mle_logistic_gradient <- function(beta, x_a, x_b, d, p_mat) {
  part1 <- colSums(x_a)
  pr_b <- as.numeric(plogis(x_b %*% beta))
  part2 <- colSums(d*x_b*pr_b)
  return(part1-part2)
}

gee_logistic <- function(beta, x_a, x_b, d, p_mat) {
  pr_a <- as.numeric(plogis(x_a %*% beta))
  res <- colSums(x_a/pr_a) - colSums(x_b*d)
  return(res)
}

# Generate some example data
set.seed(123)
n <- 10000
x1 <- sample(x = c("A", "B", "C"), size =n, replace = T, prob = c(0.70, 0.20, 0.10))
x2 <- sample(x = c("N", "A"), size = n, replace = T, prob = c(0.60, 0.40))
p <- plogis(1 - (x1=="A") + (x2 == "A"))
r <- rbinom(n, size = 1, prob = p)
y <- 2 + 1.5*(x1 == "B") - 2*(x1 == "C") + 3*(x2 == "N") + rnorm(n)
df_pop <- data.frame(x1,x2,r,y)

sample_a <- subset(df_pop, r == 1)
sample_b <- df_pop[sample(1:nrow(df_pop), 1000), ]
sample_b$d <- n/1000
sample_b_svy <- svydesign(ids=~1, data = sample_b, weights = ~ d)

res <- nonprob(selection = ~x1+x2,
               data = sample_a,
               svydesign = sample_b_svy,
               target = ~y)
summary(res)

mle_logistic(beta = c(1,1,1),
             x_a = model.matrix(~x1+x2,data=sample_a),
             x_b = model.matrix(~x1+x2,data=sample_b),
             d = sample_b$d)

mle_logistic_gradient(beta = c(1,0.5,0.3),
                      x_a = model.matrix(~x1+x2,data=sample_a),
                      x_b = model.matrix(~x1+x2,data=sample_b),
                      d = sample_b$d)

opt_result <- optim(par = c(1,1,1),
                    fn = mle_logistic,
                    gr = mle_logistic_gradient,
                    x_a = model.matrix(~x1+x2,data=sample_a),
                    x_b = model.matrix(~x1+x2,data=sample_b),
                    d = sample_b$d)
opt_result$par

gee_logistic(opt_result$par,
             x_a = model.matrix(~x1+x2,data=sample_a),
             x_b = model.matrix(~x1+x2,data=sample_b),
             d = sample_b$d)

multiroot(f = gee_logistic,
          start = opt_result$par,
          x_a = model.matrix(~x1+x2,data=sample_a),
          x_b = model.matrix(~x1+x2,data=sample_b),
          d = sample_b$d)$root

res <- nonprob(selection = ~x1+x2,
               data = sample_a,
               svydesign = sample_b_svy,
               target = ~target,
               control_selection = controlSel(est_method_sel = "gee"))

res$selection$coefficients
