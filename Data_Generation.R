# Set working directory ---------------------------------------------------
rm(list=ls())
library(tidyverse)
library(nleqslv)
setwd("F:/yifei/UCL stats for clinical trial/Simulation with Rafael/Simulations")

# get a list of all files in functions folder
pathnames <- list.files(pattern="[.]R$", path="function", full.names=TRUE)
# source all of the files in functions folder
sapply(pathnames, FUN=source)

# Modify these parameters to generate different data
control_median = 6
control_variance = 16
control_zero_prop = 0.14
intervention_median = 4
intervention_variance = 16
intervention_zero_prop = 0.2

control_param <- get_parameters(target_median = control_median, 
                                target_variance = control_variance, 
                                zero_prop=control_zero_prop)
intervention_param <- get_parameters(target_median = intervention_median,
                                     target_variance = intervention_variance,
                                     zero_prop=intervention_zero_prop)
# create a test sample of 20000 to confirm distribution
sample_size <- 20000
df <- data.frame(treatment = rep(NA, sample_size))
df <- df %>% 
  mutate(treatment=rbinom(sample_size, size=1, prob=0.5))
df <- df %>% 
  mutate(outcome = ifelse(treatment==0,
                          generate_outcome_data(sample_size, control_param[1], control_param[2], control_zero_prop),
                          generate_outcome_data(sample_size, intervention_param[1], intervention_param[2], intervention_zero_prop)))

write_csv(df, "data/sample_data.csv")

