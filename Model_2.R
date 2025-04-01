# Title
# Sub-title
# DD/MM/YYYY
# Purpose:

# Set working directory ---------------------------------------------------
rm(list=ls())
library(tidyverse)
library(nleqslv)
library(glue)
library(rsimsum)
library(glue)
setwd("F:/yifei/UCL stats for clinical trial/Simulation with Rafael/Simulations")

# get a list of all files in functions folder
pathnames <- list.files(pattern="[.]R$", path="function", full.names=TRUE)
# source all of the files in functions folder
sapply(pathnames, FUN=source)

# Design parameters of the model
control_median = 6
control_variance = 16
control_zero_prop = 0.14
intervention_median = 2
intervention_variance = 16
intervention_zero_prop = 0.24


#TEMPORARY 
mean_intervention=0
mean_control = 0


control_param <- get_parameters(target_median = control_median, 
                                target_variance = control_variance, 
                                zero_prop=control_zero_prop)
intervention_param <- get_parameters(target_median = intervention_median,
                                     target_variance = intervention_variance,
                                     zero_prop=intervention_zero_prop)



# keep record of power for each sample size
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
num_iter <- 1000
for (sample_size in seq(20, 160, 5)){
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
  res <- data.frame(
    num = numeric(),
    method = character(),
    estimand = character(),
    b = numeric(),
    se = numeric(),
    true_b = numeric(),
    p_value = numeric()
  )
  
  for (i in 1:num_iter){
    # data generation
    df <- data.frame(treatment = rep(NA, sample_size))
    df <- df %>% 
      mutate(treatment=rbinom(sample_size, size=1, prob=0.5))
    df <- df %>% 
      mutate(outcome = ifelse(treatment==0,
                              generate_outcome_data(sample_size, control_param[1], control_param[2], control_zero_prop),
                              generate_outcome_data(sample_size, intervention_param[1], intervention_param[2], intervention_zero_prop)))
    
    # using linear regression
    model <- lm(outcome ~ treatment, data=df)
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res <- rbind(res, data.frame(
      num = i,
      method = "LR",
      estimand = "mean difference",
      b = point_estimate,
      se = se,
      true_b = mean_intervention-mean_control,
      p_value = p_value
    ))
    
    # using linear regression on log-transformed data
    df$log_outcome <- log(df$outcome)
    model <- lm(log_outcome ~ treatment, data=df)
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res <- rbind(res, data.frame(
      num = i,
      method = "Log_LR",
      estimand = "log median ratio",
      b = point_estimate,
      se = se,
      true_b = log(intervention_median/control_median),
      p_value = p_value
    ))
    
    # using Gaussian glm with log link
    model <- glm(outcome ~ treatment, data=df, family=gaussian(link="log"))
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res <- rbind(res, data.frame(
      num = i,
      method = "GLM_Gaussian",
      estimand = "log mean ratio",
      b = point_estimate,
      se = se,
      true_b = log(mean_intervention/mean_control),
      p_value = p_value
    ))
    
    # using Gamma glm with log link
    model <- glm(outcome ~ treatment, data=df, family=Gamma(link="log"))
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res <- rbind(res, data.frame(
      num = i,
      method = "GLM_Gamma",
      estimand = "log mean ratio",
      b = point_estimate,
      se = se,
      true_b = log(mean_intervention/mean_control),
      p_value = p_value
    ))
    
    # using Ordinal Logistic Regression
    dd <- datadist(df)
    options(datadist = "dd")
    model <- orm(outcome ~ treatment, data=df, se.fit=TRUE)
    summ <- summary(model)
    
    point_estimate = summ['treatment','Effect']
    se = summ['treatment', "S.E."]
    z_score <- point_estimate / se
    # manually compute p-value since it's not stored
    # checked it gives the same as output
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    res <- rbind(res, data.frame(
      num = i,
      method = "OLR",
      estimand = "log odds?",
      b = point_estimate,
      se = se,
      true_b = 0,
      p_value = p_value
    ))
  }
  
  # use multisimsum to get the power and MC error for power
  s1 <-  multisimsum(data=res, par="method", true="true_b", estvarname="b", se="se")
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

write_csv(result, glue("F:/yifei/UCL stats for clinical trial/Simulation with Rafael/Simulations/result/model_2_iterations_{num_iter}.csv"))


