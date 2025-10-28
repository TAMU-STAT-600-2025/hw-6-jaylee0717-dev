
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library("microbenchmark")

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
# Test 1
a1 <- 3; lambda1 <- 1
r_soft1 <- soft(a1, lambda1)
cpp_soft1 <- soft_c(a1, lambda1)
cat("Test 1:", all.equal(r_soft1, cpp_soft1), "\n")

# Test 2
a2 <- -0.5; lambda2 <- 1
r_soft2 <- soft(a2, lambda2)
cpp_soft2 <- soft_c(a2, lambda2)
cat("Test 2:", all.equal(r_soft2, cpp_soft2), "\n")

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
set.seed(42)

# Test 1
X <- matrix(rnorm(20), nrow = 5)
Y <- rnorm(5)
beta <- rnorm(4)
lambda <- 0.1

r_lasso1 <- lasso(X, Y, beta, lambda)
cpp_lasso1 <- lasso_c(X, Y, beta, lambda)
cat("Test 1:", all.equal(as.numeric(r_lasso1), cpp_lasso1), "\n")

# Test 2 (different random inputs)
X2 <- matrix(rnorm(18), nrow = 6)
Y2 <- rnorm(6)
beta2 <- runif(3)
lambda2 <- 0.5

r_lasso2 <- lasso(X2, Y2, beta2, lambda2)
cpp_lasso2 <- lasso_c(X2, Y2, beta2, lambda2)
cat("Test 2:", all.equal(as.numeric(r_lasso2), cpp_lasso2), "\n")

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################
set.seed(99)
X3 <- matrix(rnorm(30), nrow = 10, ncol = 3)
Y3 <- rnorm(10)
lambda3 <- 0.1

# Test 1: zero start
fit_r1 <- fitLASSOstandardized(X3, Y3, lambda3)
fit_c1 <- fitLASSOstandardized_c(X3, Y3, lambda3, rep(0, ncol(X3)))
cat("Test 1:", all.equal(as.numeric(fit_r1$beta), as.numeric(fit_c1)), "\n")

# Test 2: random start
beta_start <- runif(ncol(X3))
fit_r2 <- fitLASSOstandardized(X3, Y3, lambda3, beta_start)
fit_c2 <- fitLASSOstandardized_c(X3, Y3, lambda3, beta_start)
cat("Test 2:", all.equal(as.numeric(fit_r2$beta), as.numeric(fit_c2)), "\n")



# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

# Benchmark
bench1 <- microbenchmark(
  R = fitLASSOstandardized(X2, Y2, lambda3),
  Cpp = fitLASSOstandardized_c(X2, Y2, lambda3, rep(0, ncol(X2))),
  times = 10
)
print(bench1)

bench2 <- microbenchmark(
  R = fitLASSOstandardized(X3, Y3, lambda3),
  Cpp = fitLASSOstandardized_c(X3, Y3, lambda3, rep(0, ncol(X3))),
  times = 10
)
print(bench2)



# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################
set.seed(99)
X4 <- matrix(rnorm(100), nrow = 20, ncol = 5)
Y4 <- rnorm(20)
lambda_seq <- exp(seq(log(1), log(0.01), length.out = 5))

# Test 1
fit_seq_R1 <- fitLASSOstandardized_seq(X4, Y4, lambda_seq)
fit_seq_C1 <- fitLASSOstandardized_seq_c(X4, Y4, lambda_seq)
cat("Test 1: beta matrix equality =", 
    all.equal(as.numeric(fit_seq_R1$beta_mat), as.numeric(fit_seq_C1)), "\n")

# Test 2: smaller problem, different random data
set.seed(101)
X5 <- matrix(rnorm(60), nrow = 12, ncol = 5)
Y5 <- rnorm(12)
lambda_seq2 <- exp(seq(log(0.5), log(0.05), length.out = 4))

fit_seq_R2 <- fitLASSOstandardized_seq(X5, Y5, lambda_seq2)
fit_seq_C2 <- fitLASSOstandardized_seq_c(X5, Y5, lambda_seq2)
cat("Test 2: beta matrix equality =", 
    all.equal(as.numeric(fit_seq_R2$beta_mat), as.numeric(fit_seq_C2)), "\n")



# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################
bench3 <- microbenchmark(
  R = fitLASSOstandardized_seq(X4, Y4, lambda_seq),
  Cpp = fitLASSOstandardized_seq_c(X4, Y4, lambda_seq),
  times = 5
)
print(bench3)

bench4 <- microbenchmark(
  R = fitLASSOstandardized_seq(X5, Y5, lambda_seq2),
  Cpp = fitLASSOstandardized_seq_c(X5, Y5, lambda_seq2),
  times = 5
)
print(bench4)


# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
