# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  Y <- as.vector(Y)
  n <- nrow(X)
  p <- ncol(X)
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  # [ToDo] Center and scale X
  Xmeans <- colMeans(X)
  Xcentered <- X - matrix(Xmeans, n, p, byrow = TRUE)
  weights <- sqrt(colSums(Xcentered^2 / n))
  Xtilde <- sweep(Xcentered, 2, weights, FUN = "/")
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}


# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  return(sign(a) * max(abs(a) - lambda, 0))
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  n <- nrow(Xtilde)
  return(crossprod(Ytilde - Xtilde %*% beta) / (2*n) + lambda * sum(abs(beta)))
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  Ytilde <- as.vector(Ytilde)
  if (nrow(Xtilde) != length(Ytilde)){
    stop("Number of rows in Xtilde and Ytilde must match")
  }
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  
  #[ToDo]  Check that lambda is non-negative
  if (lambda < 0){
    stop("Lambda should be nonnegative.")
  }
  
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if (is.null(beta_start)){
    beta_start <- rep(0, p)
  }
  else{
    if (length(beta_start) != p){
      stop("Length of beta_start must match p.")
    }
  }
  
  
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  fval_old <- lasso(Xtilde = Xtilde, Ytilde = Ytilde, beta = beta_start, lambda = lambda)
  beta <- beta_start
  beta_new <- beta_start
  resid <- Ytilde - Xtilde %*% beta
  
  # Convergence check at the tail
  repeat{
    # Update Beta by coordinates
    for (j in 1:p) {
      # Update Beta
      partial_resid <- resid + (Xtilde[, j] * beta[j])
      beta_new[j] <- soft((1 / n) * crossprod(Xtilde[, j], partial_resid), lambda)
      
      # Recalculate Residuals and store new beta into old
      resid <- resid - Xtilde[, j] * (beta_new[j] - beta[j])
      beta <- beta_new
    }

    
    fval_new <- lasso(Xtilde, Ytilde, beta_new, lambda)
    # Stopping condition
    if (abs(fval_old - fval_new) < eps) {
      break
    }
    
    # Update the values for the next iteration
    fval_old <- fval_new
  }
  fmin <- fval_new
  
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}


# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  Ytilde <- as.vector(Ytilde)
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != length(Ytilde)){
    stop("Number of rows in Xtilde and Ytilde must match")
  }
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.

  if (!is.null(lambda_seq)) {
    # Filter and sort
    lambda_seq <- lambda_seq[lambda_seq >= 0]
    lambda_seq <- sort(lambda_seq, decreasing = TRUE)
    
    # Check for emptiness
    if (length(lambda_seq) == 0) {
      warning("None of the supplied lambda values were non-negative. A default sequence will be generated.")
      # Set lambda_seq back to NULL to trigger the default calculation block.
      lambda_seq <- NULL
    }
  }
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  if (is.null(lambda_seq)){
    lambda_max <- max(abs(crossprod(Xtilde, Ytilde) / n))
    lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  n_lambda <- length(lambda_seq)
  beta_mat <- matrix(0, ncol = n_lambda, nrow = p)
  fmin_vec <- rep(0, n_lambda)
  
  # Main Loop
  beta_start <- rep(0, p)
  for (ind in 1:n_lambda){
    fit <- fitLASSOstandardized(Xtilde = Xtilde, Ytilde = Ytilde, lambda = lambda_seq[ind],
                                          beta_start = beta_start, eps = eps)
    fmin_vec[ind] <- fit$fmin
    beta_mat[, ind] <- fit$beta
    beta_start <- fit$beta

  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  std_res <- standardizeXY(X, Y)
  Xtilde <- std_res$Xtilde
  Ytilde <- std_res$Ytilde
  Ymean <- std_res$Ymean
  Xmeans <- std_res$Xmeans
  weights <- std_res$weights
  
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  fit <- fitLASSOstandardized_seq(Xtilde = Xtilde, Ytilde = Ytilde, lambda_seq = lambda_seq, 
                          n_lambda = n_lambda, eps = eps)
  
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  lambda_seq_used <- fit$lambda_seq
  beta_mat <- fit$beta_mat / weights
  beta0_vec <- as.vector(Ymean - crossprod(Xmeans, beta_mat))
  
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq_used, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  n <- nrow(X)
  
  # [ToDo] Fit Lasso on original data using fitLASSO
  original_fit <- fitLASSO(X = X, Y = Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  lambda_seq <- original_fit$lambda_seq
  n_lambda <- length(lambda_seq)
  
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if (is.null(fold_ids)){
    fold_ids <- sample(rep(1:k, length.out = n), size = n)
  }
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  unique_folds <- sort(unique(fold_ids))
  k_actual <- length(unique_folds)
  fold_errors <- matrix(NA, nrow = k_actual, ncol = n_lambda)


  for (i in 1:k_actual){
    # Identify training and validation data
    j <- unique_folds[i]
    val_indices <- which(fold_ids == j)
    train_indices <- which(fold_ids != j)
    
    X_train <- X[train_indices, ]
    Y_train <- Y[train_indices]
    X_val <- X[val_indices, ]
    Y_val <- Y[val_indices]
    # Fit Lasso on training data, make predictions
    fold_fit <- fitLASSO(X = X_train, Y = Y_train, lambda_seq = lambda_seq, eps = eps)
    
    # Make predictions on the validation set for the entire path
    # beta0_vec needs to be replicated for each observation in the validation set
    beta0_mat <- matrix(fold_fit$beta0_vec, nrow = length(val_indices), ncol = n_lambda, byrow = TRUE)
    predictions <- X_val %*% fold_fit$beta_mat + beta0_mat
    
    # Compute fold errors and store
    fold_errors[i, ] <- colMeans((Y_val - predictions)^2)
  }
  # Compute CVm and cvse and store
  cvm <- colMeans(fold_errors)
  cvse <- apply(fold_errors, 2, sd) / sqrt(k_actual)
  
  
  # [ToDo] Find lambda_min
  min_cvm_index <- which.min(cvm)
  lambda_min <- lambda_seq[min_cvm_index]
  
  # [ToDo] Find lambda_1SE
  one_se_threshold <-cvm[min_cvm_index] + cvse[min_cvm_index]
  lambda_1se <- max(lambda_seq[cvm <= one_se_threshold])
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = original_fit$lambda_seq, beta_mat = original_fit$beta_mat, beta0_vec = original_fit$beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

