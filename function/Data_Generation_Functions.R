get_log_normal_parameters <- function(median, variance) {
  # Purpose: Return the parameters (meanlog and sdlog) of the log-normal
  #          distribution given the median and the variance
  
  # Arguments: median - the median of the log-normal distribution
  #            variance - the variance of the log-normal distribution
  
  # Returns: vector of length 2 containing the meanlog and sdlog
  
  meanlog = log(median)
  sdlog = sqrt(log((1+sqrt(1+4*variance/median^2))/2))
  return(c(meanlog, sdlog))
}

equations <- function(params, zero_prop=0.0, target_median=0, target_variance=1) {
  # Purpose: Return the equations to solve for parameters of inflated log-normal
  #          distribution
  
  mu <- params[1]
  sigma <- params[2]
  
  # equation for the median
  composite_median = qlnorm((0.5-zero_prop)/(1-zero_prop), meanlog=mu, sdlog=sigma)
  
  # equation for variance
  composite_variance = (1-zero_prop)*(exp(sigma^2)-1+zero_prop)*exp(2*mu+sigma^2)
  
  c(
    composite_median-target_median,
    composite_variance-target_variance
  )
}

get_parameters <- function(target_median=0, target_variance=1, zero_prop=0.0) {
  # Purpose: Generate the parameters for the log-normal distribution of the given model.
  #
  # Arguments: target_median - the target median of the distribution
  #            target_variance - the target variance of the distribution
  #            zero_prop - the proportions of zeros
  #            if zero_prop is 0, then it uses direct computations of log-normal parameters,
  #            otherwise, it numerically solves equations to find the meanlog and sdlog
  
  if (zero_prop==0) {
    params <- get_log_normal_parameters(target_median, target_variance)
  } else {
    initial_guess = get_log_normal_parameters(target_median, target_variance)
    solution <- nleqslv(initial_guess,
                        function(params) equations(params,
                                                   zero_prop=zero_prop,
                                                   target_median=target_median,
                                                   target_variance=target_variance))
    params <- c(solution$x[1], solution$x[2])
  }
  params
}

generate_outcome_data <- function(n, meanlog, sdlog, zero_prop=0.0) {
  # Purpose: Generate a sample of size n from the zero_inflated log-normal model
  #          For analysis simplicity, 
  #
  # Arguments: n - the number of samples to generate
  #            meanlog - the meanlog parameter for the log-normal
  #            sdlog - the sdlog parameter for the log-normal
  #            zero_prop - the proportions of zeros
  
  ifelse(runif(n)<zero_prop, 0.01, rlnorm(n, meanlog, sdlog))
}