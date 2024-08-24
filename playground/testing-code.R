# Create a larger sample dataset
library(nonprobsvy)
library(nonprobsvyMC)
library(data.table)
library(ggplot2)

set.seed(123)
data <- data.frame(
  education = factor(sample(c("High School", "Bachelor's", "Master's", "PhD"), N, replace = TRUE)),
  income = factor(sample(c("Low", "Medium", "High"), N, replace = TRUE, prob = c(0.3, 0.5, 0.2))),
  age = factor(cut(rnorm(N, mean = 40, sd = 15),
                   breaks = c(-Inf, 25, 35, 45, 55, Inf),
                   labels = c("18-25", "26-35", "36-45", "46-55", "55+")))
)

# Define classification matrices with entries at least 0.1
class_mat <- list(
  education = matrix(
    c(0.7, 0.1, 0.1, 0.1,
      0.1, 0.6, 0.2, 0.1,
      0.1, 0.2, 0.5, 0.2,
      0.1, 0.1, 0.2, 0.6),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("High School", "Bachelor's", "Master's", "PhD"),
                    c("High School", "Bachelor's", "Master's", "PhD"))
  ),
  income = matrix(
    c(0.7, 0.2, 0.1,
      0.2, 0.6, 0.2,
      0.1, 0.3, 0.6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("Low", "Medium", "High"), c("Low", "Medium", "High"))
  ),
  age = matrix(
     c(0.5, 0.2, 0.1, 0.1, 0.1,
       0.2, 0.4, 0.2, 0.1, 0.1,
       0.1, 0.2, 0.4, 0.2, 0.1,
       0.1, 0.1, 0.2, 0.4, 0.2,
       0.1, 0.1, 0.1, 0.2, 0.5),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(c("18-25", "26-35", "36-45", "46-55", "55+"),
                    c("18-25", "26-35", "36-45", "46-55", "55+"))
  )
)


## selection to non-probability sample
set.seed(2024)
misclassified_data <- misclassify(data, class_mat)
R <- with(data, plogis(-2 -2*(income == "High") - 2*(age == "46-55") + (education == "PhD")))
Y <- with(data, plogis(1 + 2*(income == "High") - 2*(age == "26-35") -(education == "PhD")))
data$Y <- rbinom(N, 1, Y)

x_class_mat <- class_mat
x_class_mat$education <- prop.table(table(data$education, misclassified_data$education), margin = 2)
x_class_mat$income <- prop.table(table(data$income, misclassified_data$income), margin = 2)
x_class_mat$age <- prop.table(table(data$age, misclassified_data$age), margin = 2)

results_sim <- list()

for (b in 1:50) {
  set.seed(b)
  print(b)
  data$R <- rbinom(N, 1, R)
  sample_a <- subset(data, R == 1)
  sample_b <- data[sample(1:N, 1000), ]
  sample_b$d <- N/1000
  sample_b <- svydesign(ids=~1,data=sample_b, weights=~d)

  res_ipw <- nonprob(selection = ~ income + age + education,
                     target= ~ Y,
                     data = sample_a,
                     svydesign = sample_b)

  res_mi <- nonprob(outcome = Y ~ income + age + education,
                    data = sample_a,
                    svydesign = sample_b,
                    family_outcome = "binomial")

  ## with misclassified
  misclassified_data$R <- data$R
  misclassified_data$Y <- data$Y

  mis_sample_a <- subset(misclassified_data, R == 1)

  mis_res_ipw <- nonprob(selection = ~ income + age + education,
                         target= ~ Y,
                         data = mis_sample_a,
                         svydesign = sample_b)

  mis_res_mi <- nonprob(outcome = Y ~ income + age + education,
                        data = mis_sample_a,
                        svydesign = sample_b,
                        family_outcome = "binomial")

  mis_simex <- correct_misclass(obj = mis_res_mi,
                   x_class_mat = x_class_mat,
                   method = "simex",
                   svy = sample_b
                   , control = list(simex.lambda = c(0.05, 0.1, 0.25, 0.5, 1, 1.5), simex.boot = 50)
  )
  names(mis_simex) <- c("type", "mean", "SE", "lower_bound", "upper_bound")

  ## imputation based on stratified sample

  results_sim[[b]] <- rbind(
    data.frame(type = "naive", mean = mean(sample_a$Y), SE=NA, lower_bound=NA, upper_bound=NA),
    data.frame(type="ipw true", res_ipw$output, res_ipw$confidence_interval),
    data.frame(type="mi true", res_mi$output, res_mi$confidence_interval),
    data.frame(type="ipw mis", mis_res_ipw$output, mis_res_ipw$confidence_interval),
    data.frame(type="mi mis", mis_res_mi$output, mis_res_mi$confidence_interval),
    mis_simex
  )
}

results_sim_df <- rbindlist(results_sim, idcol = "b")
results_sim_df[, type:=factor(type, c("naive", "mi true", "ipw true",
                                      "mi mis", "ipw mis", "Y_lin", "Y_log", "Y_qua"))]

ggplot(data = results_sim_df, aes(x = type, y=mean)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.1) +
  geom_hline(yintercept = mean(data$Y), col = "red")


