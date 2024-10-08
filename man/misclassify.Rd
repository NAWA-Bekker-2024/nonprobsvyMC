% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misclassify.R
\name{misclassify}
\alias{misclassify}
\title{Accounting for Misclassification in nonprobsvy Auxiliary Variables}
\usage{
misclassify(data, class_mat)
}
\arguments{
\item{data}{A data frame containing the variables to be misclassified.
The variables should be factors or character vectors.}

\item{class_mat}{A named list of classification matrices.
The names should correspond to column names in the data frame.
Each matrix should be square, with dimensions equal to the number of levels
in the corresponding variable. The rows and columns should be named with the
levels of the corresponding variable. The values represent probabilities of classification.}
}
\value{
A data frame with the same structure as the input, but with
specified variables misclassified according to the provided matrices.
}
\description{
This function extends the nonprobsvy package by allowing misclassification
in auxiliary variables. It takes a data frame and a list of classification
matrices, and returns a new data frame with specified variables misclassified
according to the provided matrices.
}
\details{
The function performs the following steps:
\enumerate{
\item Checks if the input is valid (data frame and list of matrices).
\item Verifies that all specified columns exist in the data frame.
\item For each specified column, it checks if the classification matrix
matches the levels of the variable and adjusts if necessary.
\item Applies the misclassification using the corresponding classification matrix.
\item Returns the modified data frame with misclassified variables.
}

The classification matrices are normalized so that each row sums to 1.
If a classification matrix contains levels not present in the data,
those rows and columns are removed. If the data contains levels not
present in the classification matrix, those levels are added with
a diagonal probability of 1 (no misclassification for those levels).
}
\examples{
# Create a larger sample dataset
set.seed(123)
data <- data.frame(
  Education = factor(sample(c("High School", "Bachelor's", "Master's", "PhD"), 1000, replace = TRUE)),
  Income = factor(sample(c("Low", "Medium", "High"), 1000, replace = TRUE, prob = c(0.3, 0.5, 0.2))),
  Age = factor(cut(rnorm(1000, mean = 40, sd = 15),
                   breaks = c(-Inf, 25, 35, 45, 55, Inf),
                   labels = c("18-25", "26-35", "36-45", "46-55", "55+")))
)

# Define classification matrices with entries at least 0.1
class_mat <- list(
  Education = matrix(
    c(0.7, 0.1, 0.1, 0.1,
      0.1, 0.6, 0.2, 0.1,
      0.1, 0.2, 0.5, 0.2,
      0.1, 0.1, 0.2, 0.6),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("High School", "Bachelor's", "Master's", "PhD"),
                    c("High School", "Bachelor's", "Master's", "PhD"))
  ),
  Income = matrix(
    c(0.7, 0.2, 0.1,
      0.2, 0.6, 0.2,
      0.1, 0.3, 0.6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("Low", "Medium", "High"), c("Low", "Medium", "High"))
  ),
  Age = matrix(
    c(0.6, 0.2, 0.1, 0.1, 0.0,
      0.2, 0.5, 0.2, 0.1, 0.0,
      0.1, 0.2, 0.4, 0.2, 0.1,
      0.0, 0.1, 0.2, 0.5, 0.2,
      0.0, 0.0, 0.1, 0.3, 0.6),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(c("18-25", "26-35", "36-45", "46-55", "55+"),
                    c("18-25", "26-35", "36-45", "46-55", "55+"))
  )
)

# Apply misclassification
misclassified_data <- misclassify(data, class_mat)

# Compare original and misclassified data for Education
education_comparison <- table(Original = data$Education, Misclassified = misclassified_data$Education)
print("Education Misclassification Comparison:")
print(education_comparison)

# Calculate misclassification rates for Education
education_rates <- prop.table(education_comparison, margin = 1)
print("\nEducation Misclassification Rates:")
print(round(education_rates, 2))

# Compare with the original classification matrix for Education
print("\nOriginal Education Classification Matrix:")
print(round(class_mat$Education, 2))

# Repeat for Income and Age if desired

}
