# Computing the integral term in equation 3.11 ---------------------------------------
# Compute Fy(dx) in Equation (3.9)
fisher.esd = function (y1, y2) 
{c2 = y2
c1 = y1
h = sqrt(c1 + c2 - c1 * c2)
a = ((1 - h)^2)/(1 - c2)^2
b = ((1 + h)^2)/(1 - c2)^2
function(x) {
  result = rep(0, length(x))
  result[x > a & x < b] = ((1 - c2) * sqrt((b - x) * (x - a)))/(2 * pi * x * (c1 + x * c2))
  return(result)
}
}

# Computes the support for the Fisher ESD
construct.fisher.support = function (y1, y2) 
{ c2 = y2
c1 = y1
h = sqrt(c1 + c2 - c1 * c2)
a = ((1 - h)^2)/(1 - c2)^2
b = ((1 + h)^2)/(1 - c2)^2
return(c(a, b))
}

# Calculates f*(x)
function.prod = function(f1,f2){
  return(function(x){ return(f1(x)*f2(x))})
}  

### Calculates the second term after minus in the test statistic T tilda
calculate.expected.trace = function (y1, y2) 
{ asymptotic.pdf = fisher.esd(y1, y2)
integrand = function.prod(function(x) {
  (1 - x)^2 + (1 - 1/x)^2
}, asymptotic.pdf)
asymptotic.supports = construct.fisher.support(y1, y2)
safe_integral = safely(integrate)
integral = exec(safe_integral, integrand, `!!!`(asymptotic.supports))[[1]]
return(integral)
}
# Calculating the T(S1,S2) in (3.11) --------------------------------------
# Covariance distance estimator between Sigma1 and Sigma2 for tau
covariance.distance.estimator <- function (data, epsilon, f) # epsilon is a small value used for stability
{                                                   # f is the distance function between two covariances
  products = purrr::map(as.data.frame(t(data)), ~.x %*% t(.x)) # computes the product of each row with itself
  forward.cumsum = purrr::accumulate(products, ~.x + .y) # cumulative sum of the products in the previous step
  backward.cumsum = purrr::map(forward.cumsum, ~-.x + forward.cumsum[[nrow(data)]]) # cumulative sum in reverse order
  function(tau) { # calculates the covariance distance at a specific time point tau
    sigma1 = (1/tau) * forward.cumsum[[tau]] + epsilon * diag(ncol(data)) # Sigma_1   
    sigma2 = (1/(nrow(data) - tau)) * backward.cumsum[[tau]] + epsilon * diag(ncol(data)) # Sigma_2
    output = tryCatch(cov.dist(sigma1, sigma2, f), error = function(error_message) { # calculates the cov_distance based on a function f
      return(NA)# return NA if error
    })
    return(output)
  }
}

# Computes the distance between two covariance matrices using eigenvalues - equation 2.5
cov.dist= function(sigma1, sigma2, f){
  A = geigen::geigen(sigma2, sigma1, symmetric=TRUE)$values 
  return(rlang::exec(f,(A-1)^2) + rlang::exec(f,(1/A-1)^2))
}


# Calculating mu(gamma) and sigma_square ----------------------------------
# Calculates mu(gamma)
asymptotic.bias = function (y1, y2) 
{ h = sqrt(y1 + y2 - y1 * y2)
K_1 = 2 * h * (1 + h^2)/(1 - y2)^4 - 2 * h/(1 - y2)^2
J_1 = 2 * h * (1 + h^2)/(1 - y1)^4 - 2 * h/(1 - y1)^2
return(2 * (h^2 - y2^2)/(1 - y2)^4 + 2 * K_1 * y2/h + 2 * 
         (h^2 - y1^2)/(1 - y1)^4 + 2 * J_1 * y1/h)
}

# calculate the sigma square function, Wrong formula in the return
asymptotic.variance = function (y1, y2) 
{ h = sqrt(y1 + y2 - y1 * y2)
K_21 = 2 * h * (1 + h^2)/(1 - y2)^4 - 2 * h/(1 - y2)^2
K_22 = 2 * h * (1 + h^2)/(1 - y1)^4 - 2 * h/(1 - y1)^2
K_31 = h^2/(1 - y2)^4
K_32 = h^2/(1 - y1)^4
J_1 = -2 * (1 - y2)^2
J_2 = (1 - y2)^4
var_x = K_21^2 + 2 * K_31^2
var_y = K_22^2 + 2 * K_32^2
cov_xy = J_1 * K_21/h + J_1 * K_21/(h * (h^2 - 1)) + (-J_1 * 
                                                        K_31 * (h^2 + 1)/h^2) + (-J_1 * K_31/(h^2 * (h^2 - 1))) 
J_2 * K_21 * h/(h^2 - 1)^3 + J_2 * K_31/h^2 + J_2 * + # after K21 got rid of 2
  K_31 * ((1 - 3 * h^2)/(h^2 * (h^2 - 1)^3))
return(2 * (var_x + var_y + 1 * cov_xy)) # 2 instead of 3 and 1 time the cov instead of 2
return(limiting.var)
}

# Test Statistic ----------------------------------------------------------
# Calculates the statistic for change point detection
matrix.dist.test.stat = function (data, minseglen) 
{ p = ncol(data)
n = nrow(data)
t = seq(p + minseglen + 1, n - p - minseglen) # possible change points
estimate.covariance.distance = covariance.distance.estimator(data, 0, mean) #calculate T in eq. 3.11
test.stat = purrr::map_dbl(t, estimate.covariance.distance) # compute T for all t's
dimension.over.length = purrr::map(t, ~c(p/.x, p/(n - .x)))#calculate(??1=p/t)&(??2=p/(n-t)) for all t's
trace = map(dimension.over.length, ~exec(calculate.expected.trace, !!!.x))#for all ts calculate the integral part
values = map_lgl(trace, ~length(.x[[1]]) == 0)# check if the resultsof the integral are valid
bias = map_dbl(dimension.over.length, ~exec(asymptotic.bias, !!!.x)) # calculate ??(??) for all t's
variance = map_dbl(dimension.over.length, ~exec(asymptotic.variance, !!!.x)) # compute ??(??) for all t's
trace = map_if(trace, values, ~NA) %>% map_dbl(~.x[[1]]) # replace with NA invalid trace
bias = map_if(bias, values, ~NA)  # replace with NA invalid bias
variance = map_if(variance, values, ~NA)  # replace with NA invalid variance
test.stat = pmap_dbl(list(test.stat, trace, bias, variance), # calculate the test statistic 
                     ~(p * (..1 - ..2) - ..3)/sqrt(..4)) #
return(c(rep(NA, p + minseglen), (test.stat), rep(NA, p + minseglen)))
}
# Results -----------------------------------------------------------------
# for different possible change points
limiting.variance = function(n,p){
  t = seq(p+1,n-p)
  dimension.over.length = purrr::map(t, ~c(p/.x, p/(n-.x)))
  variance = map_dbl(dimension.over.length, ~exec(asymptotic.variance,!!!.x) ) # calcul asympt var
  return(variance)
}

# calculates from where to start the sequenc eof possible change points
calculate.minseglen = function(n, p, func, alpha=2, constraint=1/n){
  #alpha 
  if(func == "linear"){
    return(alpha*p)
  }  else if( func == "log-linear"){
    return(log(n)*p)
  } else if( func == "constrained-gradient"){
    grad = limiting.variance(n,p) %>% diff %>% abs 
    return(p + min(which(grad < constraint)))
  }
}

# Obtain a critical value
bonferoni = function(n,alpha){ qnorm(1-.05/n)}
library(tidyverse)

simulations_ryan_killick <- function(n_simulations, delta, n, p, cp, minseglen) {
  # Initialize results storage
  results <- matrix(nrow = n_simulations, ncol = 2)
  colnames(results) <- c("max_abs_result", "cp")
  
  # Main simulation loop
  for (i in 1:n_simulations) {
    # Data generation
    Y <- matrix(rnorm(n*p, 0, 1), nrow = n)
    Y[cp:n, ] <- delta * Y[cp:n, ]
    
    # Calculate test statistics with EXACT same minseglen logic
    current_minseglen <- max(4 * p, 30)
    result <- matrix.dist.test.stat(Y, current_minseglen - p)  # Exact same adjustment
    
    # Store results
    tau <- which.max(abs(result))
    results[i, ] <- c(result[tau], tau)
  }
  
  return(results)
}


# # Example usage
# set.seed(1299)
# sim_results <- simulations_ryan_killick(
#   n_simulations = 10,
#   n = 500,
#   p = 10,
#   delta = 1,
#   cp = 250,  # floor(500/2)
#   minseglen = max(4*p, 30)
# )
# 
# # View results
# head(sim_results)
