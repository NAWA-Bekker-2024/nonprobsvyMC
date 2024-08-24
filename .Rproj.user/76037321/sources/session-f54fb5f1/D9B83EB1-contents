library(tinytest)

# Test 1: Basic functionality
test_basic_functionality <- function() {
  set.seed(123)  # for reproducibility
  data <- data.frame(
    Class = factor(c("A", "B", "C", "A", "B")),
    Gender = factor(c("Male", "Female", "Male", "Female", "Male"))
  )

  class_mat <- list(
    Class = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, byrow = TRUE,
                   dimnames = list(c("A", "B", "C"), c("A", "B", "C"))),
    Gender = matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE,
                    dimnames = list(c("Male", "Female"), c("Male", "Female")))
  )

  result <- misclassify(data, class_mat)

  expect_true(all(result$Class == data$Class), "Basic functionality: Class should not change")
  expect_true(all(result$Gender == data$Gender), "Basic functionality: Gender should not change")
}

# Test 2: Complete misclassification
test_complete_misclassification <- function() {
  set.seed(123)  # for reproducibility
  data <- data.frame(
    Class = factor(c("A", "B", "C", "A", "B")),
    Gender = factor(c("Male", "Female", "Male", "Female", "Male"))
  )

  class_mat <- list(
    Class = matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE,
                   dimnames = list(c("A", "B", "C"), c("A", "B", "C"))),
    Gender = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE,
                    dimnames = list(c("Male", "Female"), c("Male", "Female")))
  )

  result <- misclassify(data, class_mat)

  expect_true(all(as.character(result$Class) != as.character(data$Class)), "Complete misclassification: All Class values should be different")
  expect_true(all(as.character(result$Gender) != as.character(data$Gender)), "Complete misclassification: All Gender values should be different")
}

# Test 3: Partial misclassification
test_partial_misclassification <- function() {
  set.seed(123)  # for reproducibility
  data <- data.frame(
    Class = factor(rep(c("A", "B", "C"), 100)),
    Gender = factor(rep(c("Male", "Female"), 150))
  )

  class_mat <- list(
    Class = matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8), nrow = 3, byrow = TRUE,
                   dimnames = list(c("A", "B", "C"), c("A", "B", "C"))),
    Gender = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE,
                    dimnames = list(c("Male", "Female"), c("Male", "Female")))
  )

  result <- misclassify(data, class_mat)

  misclassified_class <- sum(as.character(result$Class) != as.character(data$Class))
  misclassified_gender <- sum(as.character(result$Gender) != as.character(data$Gender))

  expect_true(misclassified_class > 0 && misclassified_class < nrow(data),
              "Partial misclassification: Some but not all Class values should be different")
  expect_true(misclassified_gender > 0 && misclassified_gender < nrow(data),
              "Partial misclassification: Some but not all Gender values should be different")
}

# Test 4: Error handling - Invalid input data
test_invalid_input_data <- function() {
  invalid_data <- c(1, 2, 3, 4, 5)
  class_mat <- list(matrix(c(1, 0, 0, 1), nrow = 2,
                           dimnames = list(c("A", "B"), c("A", "B"))))

  expect_error(misclassify(invalid_data, class_mat), "Input data must be a data frame")
}

# Test 5: Error handling - Invalid classification matrices
test_invalid_classification_matrices <- function() {
  data <- data.frame(Class = factor(c("A", "B", "C")))
  invalid_matrices <- list(Class = c(1, 2, 3))

  expect_error(misclassify(data, invalid_matrices), "All elements in class_mat must be matrices")
}

# Test 6: Error handling - Mismatched column names
test_mismatched_column_names <- function() {
  data <- data.frame(Class = factor(c("A", "B", "C")))
  class_mat <- list(InvalidColumn = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3,
                                           dimnames = list(c("A", "B", "C"), c("A", "B", "C"))))

  expect_error(misclassify(data, class_mat), "Some specified columns do not exist in the data frame")
}

# Test 7: Handling extra levels in classification matrix
test_extra_levels_in_matrix <- function() {
  data <- data.frame(Class = factor(c("A", "B", "A", "B")))
  class_mat <- list(
    Class = matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8), nrow = 3, byrow = TRUE,
                   dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  )

  result <- misclassify(data, class_mat)

  expect_true(all(levels(result$Class) == c("A", "B")), "Extra level 'C' should be removed from the result")
}

# Test 8: Handling missing levels in classification matrix
test_missing_levels_in_matrix <- function() {
  data <- data.frame(Class = factor(c("A", "B", "C", "A", "B")))
  class_mat <- list(
    Class = matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE,
                   dimnames = list(c("A", "B"), c("A", "B")))
  )

  result <- misclassify(data, class_mat)

  expect_true(all(levels(result$Class) == c("A", "B", "C")), "Missing level 'C' should be added to the result")
  expect_true(all(result$Class[data$Class == "C"] == "C"), "Level 'C' should not be misclassified")
}

# Run all tests
run_tests <- function() {
  test_basic_functionality()
  test_complete_misclassification()
  test_partial_misclassification()
  test_invalid_input_data()
  test_invalid_classification_matrices()
  test_mismatched_column_names()
  test_extra_levels_in_matrix()
  test_missing_levels_in_matrix()

  print("All tests completed.")
}

# Uncomment the following line to run the tests
# run_tests()
