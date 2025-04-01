# Core Likelihood Ratio Functions ------------------------------------------
compute_lr_components <- function(data_before, data_after, p, n, t) {
  # Regularized covariance matrices
  Sigma_1 <- cov(data_before) + 1e-8 * diag(p)
  Sigma_2 <- cov(data_after) + 1e-8 * diag(p)
  
  # Eigen decomposition
  F_t <- Sigma_1 %*% solve(Sigma_2)
  values <- eigen(F_t, symmetric = FALSE, only.values = TRUE)$values
  
  # Dimension ratios
  y1 <- p/(t-1)
  y2 <- p/(n-t-1)
  c1 <- (t-1)/(n-2)
  c2 <- (n-t-1)/(n-2)
  h <- sqrt(y1 + y2 - y1*y2)
  
  # Calculate u term
  u <- ((y1 + y2 - y1*y2)/(y1*y2))*log((y1 + y2)/(y1 + y2 - y1*y2)) +
    (y1*(1-y2)/(y2*(y1+y2)))*log(1-y2) +
    (y2*(1-y1)/(y1*(y1+y2)))*log(1-y1)
  
  # Calculate a and b terms
  a1 <- (1 + h^2)/((1-y2)^2)
  b1 <- 2*h/((1-y2)^2)
  a2 <- (y1*(1-y2)^2 + y2*(1+h^2))/((y1+y2)*(1-y2)^2)
  b2 <- 2*y2*h/((y1+y2)*(1-y2)^2)
  
  # Square root terms
  a2_b2_sq <- sqrt(a2^2 - b2^2)
  a1_b1_sq <- sqrt(a1^2 - b1^2)
  a2_plus_b2 <- sqrt(a2 + b2)
  a2_minus_b2 <- sqrt(a2 - b2)
  a1_plus_b1 <- sqrt(a1 + b1)
  a1_minus_b1 <- sqrt(a1 - b1)
  
  # Mu terms
  mu1 <- 0.5*log(a2_b2_sq*4*h^2/(h*(a2_plus_b2 + a2_minus_b2) - y2*(a2_plus_b2 - a2_minus_b2))^2)
  mu2 <- 0.5*log(a1_b1_sq*4*h^2/(h*(a1_plus_b1 + a1_minus_b1) - y2*(a1_plus_b1 - a1_minus_b1))^2)
  mu_n_t <- mu1 - (y2/(y1+y2))*mu2
  
  # Variance terms
  cov_term <- 2*log((a1_plus_b1 + a1_minus_b1)*(a2_plus_b2 + a2_minus_b2)/
                      ((a1_plus_b1 + a1_minus_b1)*(a2_plus_b2 + a2_minus_b2) - 
                         (a1_plus_b1 - a1_minus_b1)*(a2_plus_b2 - a2_minus_b2)))
  
  var1 <- 2*log((a2 + a2_b2_sq)/(2*a2_b2_sq))
  var2 <- 2*log((a1 + a1_b1_sq)/(2*a1_b1_sq))
  sigma_sq <- var1 + (y2/(y1+y2))^2*var2 - 2*(y2/(y1+y2))*cov_term
  
  # Calculate statistic components
  Statistic <- sum(((t-1)/2 * log(values) - ((n-2)/2 * log(c1*values + c2))))
  
  list(
    values = values,
    y1 = y1, y2 = y2, c1 = c1, c2 = c2,
    u = u, mu_n_t = mu_n_t, sigma_sq = sigma_sq,
    Statistic = Statistic
  )
}

compute_lr_statistic <- function(components, p, n) {
  with(components, {
    ((-2/(n-2))*Statistic - p*u - mu_n_t)/sqrt(abs(sigma_sq))
  })
}

# Main Analysis Function --------------------------------------------------
analyze_likelihood_ratio <- function(Y, minseglen = NULL, alpha = 0.05) {
  # Validate input
  if (!is.matrix(Y) && !is.data.frame(Y)) {
    stop("Input data must be a matrix or dataframe")
  }
  
  # Set dimensions and parameters
  p <- ncol(Y)
  n <- nrow(Y)
  if (is.null(minseglen)) minseglen <- max(4*p, 30)
  
  # Candidate change points
  t_0 <- seq(p + minseglen + 1, n - p - minseglen)
  
  # Pre-allocate results
  results <- matrix(nrow = length(t_0), ncol = 3)
  colnames(results) <- c("candidate_cp", "standardized_stat", "raw_stat")
  
  # Compute statistics for each candidate point
  for (i in seq_along(t_0)) {
    t <- t_0[i]
    components <- compute_lr_components(
      data_before = Y[1:t, , drop = FALSE],
      data_after = Y[(t+1):n, , drop = FALSE],
      p = p, n = n, t = t
    )
    
    results[i, ] <- c(
      t,
      compute_lr_statistic(components, p, n),
      components$Statistic
    )
  }
  
  # Find maximum statistic
  max_idx <- which.max(results[, "standardized_stat"])
  
  # Return comprehensive results
  list(
    test_statistics = results,
    change_point = results[max_idx, "candidate_cp"],
    max_statistic = results[max_idx, "standardized_stat"],
    critical_value = qnorm(1 - alpha/n),
    is_significant = results[max_idx, "standardized_stat"] > qnorm(1 - alpha/n),
    minseglen = minseglen,
    dimensions = c(n = n, p = p)
  )
}
plot_lr_results <- function(lr_results, start_point = NULL) {
  if (is.null(start_point)) {
    start_point <- min(lr_results$test_statistics[, "candidate_cp"])
  }
  
  # Increase right margin (bottom, left, top, right)
  par(mar = c(5, 4, 4, 5) + 0.1)  # Adjusted for better spacing
  
  plot(x = lr_results$test_statistics[, "candidate_cp"], 
       y = lr_results$test_statistics[, "standardized_stat"], 
       type = "l",
       main = "Likelihood Ratio Test Statistics", 
       xlab = "Candidate Change Point", 
       ylab = "Standardized Statistic",
       xlim = c(start_point, max(lr_results$test_statistics[, "candidate_cp"])))
  
  abline(v = lr_results$change_point, col = "red")
  abline(h = lr_results$critical_value, col = "blue", lty = 2)
  
  # Place legend outside to the right (smaller size)
  legend(x = "topright", 
         legend = c("Test Statistic", "Detected CP", "Critical Value"),
         col = c("black", "red", "blue"), 
         lty = c(1, 1, 2),
         xpd = TRUE,        # Allow plotting outside
         inset = c(-0.25, 0),  # Move right
         cex = 0.5,         # Reduce text size
         box.lty = 0)       # Remove legend box (optional)
}

# # Example Usage -----------------------------------------------------------
# 
# # Using your data
# data <- read.csv("Monthly_Prices.csv")[,4:22]

# 
# # Run analysis with minseglen = 20
# lr_results <- analyze_likelihood_ratio(data, minseglen = 20)
# 
# # Visualization starting from point 20
#plot_lr_results(lr_results, start_point = 20)
