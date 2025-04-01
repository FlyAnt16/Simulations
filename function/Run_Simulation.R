run_simulation <- function(control_median, control_variance, control_zero_prop,
                           intervention_median, intervention_variance, intervention_zero_prop,
                           num_iter, sample_size_start, sample_size_end, sample_size_step,
                           function_paths, parallel=TRUE) {
  # Purpose: Run the simulations for the given parameters
  
  
  sapply(function_paths, FUN=source)

  control_param <- get_parameters(target_median = control_median, 
                                  target_variance = control_variance, 
                                  zero_prop=control_zero_prop)
  
  intervention_param <- get_parameters(target_median = intervention_median,
                                       target_variance = intervention_variance,
                                       zero_prop=intervention_zero_prop)
  
  result <- data.frame(
    sample_size = numeric(),
    method = character(),
    p_value_power = numeric(),
    power = numeric(),
    mcse = numeric(),
    lower = numeric(),
    upper = numeric()
  )
  # start of simulation
  set.seed(1)
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  for (sample_size in seq(sample_size_start, sample_size_end, sample_size_step)){
    # keep record of results to be used for rsimsum
    # columns required:
    # 1) dataset, integer for number of samples
    # 2) method, char to indicate method
    # "LR for linear regression"
    # "Log_LR for taking log then linear regressions"
    # "GLM_Gaussian for gaussian glm with log link"
    # "OLR for ordinal logistic regression"
    # 3) b, the point estimate 
    # 4) se, the standard error of point estimate
    # 5) true_b, the true value for the point estimate
    # below are additional columns to compute power using p-value
    # 6) p-value, the p-value of the point estimate
    
    res <- foreach(i=1:num_iter, .combine=rbind,
                   .export = c("generate_outcome_data", "run_lm", "run_log_lm", "run_gaussian_log_link", "run_gamma_log_link", "run_olr"),
                   .packages=c("rms", "tidyverse")) %dopar% {
      # data generation
      df <- data.frame(treatment = rep(NA, sample_size))
      df <- df %>% 
        mutate(treatment=rbinom(sample_size, size=1, prob=0.5))
      df <- df %>% 
        mutate(outcome = ifelse(treatment==0,
                                generate_outcome_data(sample_size, control_param[1], control_param[2], control_zero_prop),
                                generate_outcome_data(sample_size, intervention_param[1], intervention_param[2], intervention_zero_prop)))
      models <- list()
      
      # using linear regression
      models[[1]] <- run_lm("outcome ~ treatment", df, i)
      
      # using linear regression on log-transformed data
      models[[2]] <- run_log_lm("outcome ~ treatment", df, i)
       
      # using Gaussian glm with log link
      models[[3]] <- run_gaussian_log_link("outcome ~ treatment", df, i)
      
      # using Gamma glm with log link
      models[[4]] <- run_gamma_log_link("outcome ~ treatment", df, i)
      
      # using Ordinal Logistic Regression
      models[[5]] <- run_olr("outcome ~ treatment", df, i)
      
      do.call(rbind, models)
                   }
    
    # use multisimsum to get the power and MC error for power
    s1 <-  multisimsum(data=res, par="method", estvarname="b", se="se")
    summ <- summary(s1)[[1]]
    summ <- summ %>% 
      filter(stat=='power') %>% 
      select(-stat) %>% 
      mutate(sample_size = sample_size) %>% 
      select(sample_size, everything()) %>% 
      rename(power=est)
    # manually compute the power using the p-value
    power_df <- res %>% 
      group_by(method) %>% 
      summarise(p_value_power = mean(p_value<0.05))
    
    # get a resultant table by mergin the previous two
    result <- rbind(result, inner_join(summ, power_df, join_by('method')))  
  }
  
  result
}

