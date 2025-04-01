rm(list=ls())
library(tidyverse)
library(nleqslv)
library(glue)
library(rsimsum)
library(rms)
library(tictoc)
setwd("F:/yifei/UCL stats for clinical trial/Simulation with Rafael/Simulations")

# get a list of all files in functions folder
pathnames <- list.files(pattern="[.]R$", path="function", full.names=TRUE)
# source all of the files in functions folder
sapply(pathnames, FUN=source)

# Design parameters of the model
control_median = 6
control_variance = 16
control_zero_prop = 0.14
intervention_median = 4
intervention_variance = 16
intervention_zero_prop = 0.2


#TEMPORARY 
mean_intervention=0
mean_control = 0


control_param <- get_parameters(target_median = control_median, 
                                target_variance = control_variance, 
                                zero_prop=control_zero_prop)
intervention_param <- get_parameters(target_median = intervention_median,
                                     target_variance = intervention_variance,
                                     zero_prop=intervention_zero_prop)

sample_size_start = 100
sample_size_end = 160
sample_size_step = 5



tic("lapply and rbind")
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
num_iter <- 500
# using lapply and rbind
for (sample_size in seq(100, 160, 5)){
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
 
  res <- lapply(seq_len(num_iter), function(i) {
    # data generation
    df <- data.frame(treatment = rep(NA, sample_size))
    df <- df %>% 
      mutate(treatment=rbinom(sample_size, size=1, prob=0.5))
    df <- df %>% 
      mutate(outcome = ifelse(treatment==0,
                              generate_outcome_data(sample_size, control_param[1], control_param[2], control_zero_prop),
                              generate_outcome_data(sample_size, intervention_param[1], intervention_param[2], intervention_zero_prop)))
    temp = data.frame(
      num = numeric(),
      method = character(),
      estimand = character(),
      b = numeric(),
      se = numeric(),
      true_b = numeric(),
      p_value = numeric()
    )
    
    # using linear regression
      model <- lm(outcome ~ treatment, data=df)
      summ <- summary(model)
      point_estimate = summ$coefficients[2,1]
      se = summ$coefficients[2,2]
      p_value <- summ$coefficients[2,4]
      temp <- rbind(temp, data.frame(
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
      temp <- rbind(temp, data.frame(
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
      temp <- rbind(temp, data.frame(
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
      temp <- rbind(temp, data.frame(
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
      temp <- rbind(temp, data.frame(
        num = i,
        method = "OLR",
        estimand = "log odds?",
        b = point_estimate,
        se = se,
        true_b = 0,
        p_value = p_value
      ))
    temp
  })
  
  res <- do.call(rbind, res)
  
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
toc()

tic("lapply and rbindlist")
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
num_iter <- 5000
# using lapply and rbindlist
for (sample_size in seq(100, 160, 5)){
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
  
  res <- lapply(seq_len(num_iter), function(i) {
    # data generation
    df <- data.frame(treatment = rep(NA, sample_size))
    df <- df %>% 
      mutate(treatment=rbinom(sample_size, size=1, prob=0.5))
    df <- df %>% 
      mutate(outcome = ifelse(treatment==0,
                              generate_outcome_data(sample_size, control_param[1], control_param[2], control_zero_prop),
                              generate_outcome_data(sample_size, intervention_param[1], intervention_param[2], intervention_zero_prop)))
    temp = data.frame(
      num = numeric(),
      method = character(),
      estimand = character(),
      b = numeric(),
      se = numeric(),
      true_b = numeric(),
      p_value = numeric()
    )
    
    # using linear regression
    model <- lm(outcome ~ treatment, data=df)
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    temp <- rbind(temp, data.frame(
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
    temp <- rbind(temp, data.frame(
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
    temp <- rbind(temp, data.frame(
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
    temp <- rbind(temp, data.frame(
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
    temp <- rbind(temp, data.frame(
      num = i,
      method = "OLR",
      estimand = "log odds?",
      b = point_estimate,
      se = se,
      true_b = 0,
      p_value = p_value
    ))
    temp
  })
  
  res <- data.table::rbindlist(res)
  
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
toc()

tic("for loops")
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
num_iter <- 5000
for (sample_size in seq(100, 160, 5)){
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
toc()

#### TODO: fixing pre-allocation, timing wasn't right
tic("pre-allocation")
# keep record of power for each sample size
result <- data.frame(
  sample_size = numeric(65),
  method = character(65),
  p_value_power = numeric(65),
  power = numeric(65),
  mcse = numeric(65),
  lower = numeric(65),
  upper = numeric(65)
)

# start of simulation
set.seed(1)
num_iter <- 5000
for (sample_size in seq(100, 160, 5)){
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
    num = numeric(25000),
    method = character(25000),
    estimand = character(25000),
    b = numeric(25000),
    se = numeric(25000),
    true_b = numeric(25000),
    p_value = numeric(25000)
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
    res[(i-1)*5+1,] <- list(
      num = i,
      method = "LR",
      estimand = "mean difference",
      b = point_estimate,
      se = se,
      true_b = mean_intervention-mean_control,
      p_value = p_value
    )
    
    # using linear regression on log-transformed data
    df$log_outcome <- log(df$outcome)
    model <- lm(log_outcome ~ treatment, data=df)
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res[(i-1)*5+2,] <- list(
      num = i,
      method = "Log_LR",
      estimand = "log median ratio",
      b = point_estimate,
      se = se,
      true_b = log(intervention_median/control_median),
      p_value = p_value
    )
    
    # using Gaussian glm with log link
    model <- glm(outcome ~ treatment, data=df, family=gaussian(link="log"))
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res[(i-1)*5+3,] <- list(
      num = i,
      method = "GLM_Gaussian",
      estimand = "log mean ratio",
      b = point_estimate,
      se = se,
      true_b = log(mean_intervention/mean_control),
      p_value = p_value
    )
    
    # using Gamma glm with log link
    model <- glm(outcome ~ treatment, data=df, family=Gamma(link="log"))
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res[(i-1)*5+4,] <- list(
      num = i,
      method = "GLM_Gamma",
      estimand = "log mean ratio",
      b = point_estimate,
      se = se,
      true_b = log(mean_intervention/mean_control),
      p_value = p_value
    )
    
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
    res[(i-1)*5+5,] <- list(
      num = i,
      method = "OLR",
      estimand = "log odds?",
      b = point_estimate,
      se = se,
      true_b = 0,
      p_value = p_value
    )
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
  result[(sample_size-100)/5*5+1 : (sample_size-100)/5*5+5, ] <- inner_join(summ, power_df, join_by('method'))
  #result <- rbind(result, inner_join(summ, power_df, join_by('method')))
}
toc()

tic("using lists")
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
num_iter <- 5000
for (sample_size in seq(100, 160, 5)){
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
  #res <- data.frame(
  #  num = numeric(25000),
  #  method = character(25000),
  #  estimand = character(25000),
  #  b = numeric(25000),
  #  se = numeric(25000),
  #  true_b = numeric(25000),
  #  p_value = numeric(25000)
  #)
  res = vector("list", length=25000)
  res_counter = 1
  
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
    
    
    res[[res_counter]] <- list(
      num = i,
      method = "LR",
      estimand = "mean difference",
      b = point_estimate,
      se = se,
      true_b = mean_intervention-mean_control,
      p_value = p_value
    )
    res_counter = res_counter + 1
    
    # using linear regression on log-transformed data
    df$log_outcome <- log(df$outcome)
    model <- lm(log_outcome ~ treatment, data=df)
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res[[res_counter]] <- list(
      num = i,
      method = "Log_LR",
      estimand = "log median ratio",
      b = point_estimate,
      se = se,
      true_b = log(intervention_median/control_median),
      p_value = p_value
    )
    res_counter = res_counter + 1
    
    # using Gaussian glm with log link
    model <- glm(outcome ~ treatment, data=df, family=gaussian(link="log"))
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res[[res_counter]] <- list(
      num = i,
      method = "GLM_Gaussian",
      estimand = "log mean ratio",
      b = point_estimate,
      se = se,
      true_b = log(mean_intervention/mean_control),
      p_value = p_value
    )
    res_counter = res_counter + 1
    
    # using Gamma glm with log link
    model <- glm(outcome ~ treatment, data=df, family=Gamma(link="log"))
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    res[[res_counter]] <- list(
      num = i,
      method = "GLM_Gamma",
      estimand = "log mean ratio",
      b = point_estimate,
      se = se,
      true_b = log(mean_intervention/mean_control),
      p_value = p_value
    )
    res_counter = res_counter + 1
    
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
    res[[res_counter]] <- list(
      num = i,
      method = "OLR",
      estimand = "log odds?",
      b = point_estimate,
      se = se,
      true_b = 0,
      p_value = p_value
    )
    res_counter = res_counter + 1
  }
  res <- data.table::rbindlist(res)
  
  
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
toc()


library(foreach)
library(doParallel)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
tic("parallel computing")

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
# using lapply and rbind
for (sample_size in seq(100, 160, 5)){
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
  
  
  res <- foreach(i=1:num_iter, .combine=rbind, .packages=c("rms", "tidyverse", "nleqslv", "rsimsum")) %dopar%{
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
    formula <- "outcome ~ treatment"
    model <- eval(bquote(lm(.(as.formula(formula)), data=df)))
    summary(model)
    
    
    model <- lm(outcome ~ treatment, data=df)
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    models[[1]] <- data.frame(
      num = i,
      method = "LR",
      estimand = "mean difference",
      b = point_estimate,
      se = se,
      true_b = mean_intervention-mean_control,
      p_value = p_value
    )
    
    # using linear regression on log-transformed data
    df$log_outcome <- log(df$outcome)
    model <- lm(log_outcome ~ treatment, data=df)
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    models[[2]] <- data.frame(
      num = i,
      method = "Log_LR",
      estimand = "log median ratio",
      b = point_estimate,
      se = se,
      true_b = log(intervention_median/control_median),
      p_value = p_value
    )
    
    # using Gaussian glm with log link
    model <- glm(outcome ~ treatment, data=df, family=gaussian(link="log"))
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    models[[3]] <- data.frame(
      num = i,
      method = "GLM_Gaussian",
      estimand = "log mean ratio",
      b = point_estimate,
      se = se,
      true_b = log(mean_intervention/mean_control),
      p_value = p_value
    )
    
    # using Gamma glm with log link
    model <- glm(outcome ~ treatment, data=df, family=Gamma(link="log"))
    summ <- summary(model)
    point_estimate = summ$coefficients[2,1]
    se = summ$coefficients[2,2]
    p_value <- summ$coefficients[2,4]
    models[[4]] <- data.frame(
      num = i,
      method = "GLM_Gamma",
      estimand = "log mean ratio",
      b = point_estimate,
      se = se,
      true_b = log(mean_intervention/mean_control),
      p_value = p_value
    )
    
    # using Ordinal Logistic Regression
    dd <- datadist(df)
    options(datadist = dd)
    model <- orm(outcome ~ treatment, data=df, se.fit=TRUE)
    summ <- summary(model)
    
    point_estimate = summ['treatment','Effect']
    se = summ['treatment', "S.E."]
    z_score <- point_estimate / se
    # manually compute p-value since it's not stored
    # checked it gives the same as output
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    models[[5]] <- data.frame(
      num = i,
      method = "OLR",
      estimand = "log odds?",
      b = point_estimate,
      se = se,
      true_b = 0,
      p_value = p_value
    )
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
toc()




mkRow <- function(nCol) {
  x <- as.list(rnorm(nCol))
  # make row mixed types by changing first column to string
  x[[1]] <- ifelse(x[[1]]>0,'pos','neg')
  names(x) <- paste('x',seq_len(nCol),sep='.')
  x
}


mkFrameList <- function(nRow,nCol) {
  d <- lapply(seq_len(nRow),function(i) {
    ri <- mkRow(nCol)
    data.frame(ri,
               stringsAsFactors=FALSE)
  })
  print(d[[1]])
  do.call(rbind,d)
  #data.table::rbindlist(d)
}


mkFrameForLoop(4,3)
y <- mkFrameList(4,3)



d <- lapply(seq_len(4),function(i) {
  ri <- mkRow(3)
  data.frame(ri,
             stringsAsFactors=FALSE)
})
print(d[[1]])
d <- do.call(rbind,d)

x <- mkRow(5)
mkFrameDataTableFor <- function(nRow, nCol) {
  v = vector("list", nRow)
  print(v)
  for (i in seq_len(nRow)) {
    v[[i]] = mkRow(nCol)
  }
  print(v)
  data.table::rbindlist(v)
}

mkFrameDataTableFor(5,2)

# store a row as a list x
# names(x) = column headers
# x[[1]], x[[2]], x[[num of cols]]
