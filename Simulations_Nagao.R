simulations_nagao <- function(n_simulations, n, p, delta, cp, minseglen) {
  # Input validation
  if(length(cp) != 1 || !is.numeric(cp)) {
    stop("cp must be a single numeric value")
  }
  
  # Initialize results matrix
  results <- matrix(nrow = n_simulations, ncol = 3)
  colnames(results) <- c("max_statistic", "cp", "detection")
  
  # Pre-compute candidate points
  minseglen <- max(4 * p, 30)
  t_0 <- seq(p + minseglen + 1, n - p - minseglen)
  
  # Main simulation loop
  for (sim in 1:n_simulations) {
    # Generate data with change point at cp
    Y <- matrix(rnorm(n * p), nrow = n)
    if(cp < n) {  # Ensure cp is within valid range
      Y[(cp+1):n, ] <- delta * Y[(cp+1):n, ]
    }
    
    # Compute test statistics safely
    test_stats <- sapply(t_0, function(t) {
      tryCatch({
        compute_nagao_test_stat(Y, t, p, n)
      }, error = function(e) NA)
    })
    
    # Store results
    if(all(is.na(test_stats))) {
      results[sim, ] <- c(NA, NA, 0)
    } else {
      max_idx <- which.max(test_stats)
      results[sim, ] <- c(
        max(test_stats, na.rm = TRUE),
        t_0[max_idx],
        as.integer(max(test_stats, na.rm = TRUE) > qnorm(1 - 0.05/n))
      )
    }
  }
  
  return(results)
}

# # Example usage
# set.seed(1299)
# n <- 500
# p <- 10
# delta <- 1
# cp <- floor(n / 2)
# minseglen <- max(4 * p, 30)
# 
# # Run simulations
# sim_results <- simulations_nagao(
#   n_simulations = 10,
#   n = n,
#   p = p,
#   delta = delta,
#   cp = cp,
#   minseglen = minseglen
# )
# 
# # View results
# head(sim_results)