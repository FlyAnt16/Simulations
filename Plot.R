# Title
# Sub-title
# DD/MM/YYYY
# Purpose:

# Set working directory ---------------------------------------------------
rm(list=ls())
library(tidyverse)
library(glue)
setwd("F:/yifei/UCL stats for clinical trial/Simulation with Rafael/Simulations")

for (model in 1:2) {
  res <- read_csv(glue("result/model_{model}_iterations_1000.csv"))
  res <- res %>% 
    mutate(p_value_mcse = sqrt(p_value_power*(1-p_value_power)/1000)) %>% 
    mutate(p_value_upper = p_value_power + qnorm(0.975)*p_value_mcse,
           p_value_lower = p_value_power - qnorm(0.975)*p_value_mcse) 
  
  plot <- ggplot(res, aes(sample_size, power, color=method)) + 
    geom_line() + 
    geom_linerange(aes(ymin=lower, ymax=upper)) +
    scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, by=0.1)) + 
    labs(title=glue("Simulation results with 1000 iterations model {model}"), x="sample size", y="power")
  ggsave(plot=plot, glue("figures/model_{model}_result_1000.jpg"))
  
  plot <- ggplot(res, aes(sample_size, p_value_power, color=method)) + 
    geom_line() + 
    geom_linerange(aes(ymin=p_value_lower, ymax=p_value_upper)) +
    scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, by=0.1)) + 
    labs(title=glue("Simulation results with 1000 iterations model {model}"), x="sample size", y="power")
  ggsave(plot=plot, glue("figures/model_{model}_pvalue_result_1000.jpg"))
}

res <- read_csv(glue("result/significance.csv"))
plot <- ggplot(res, aes(sample_size, p_value_power, color=method)) + 
  geom_line() + 
  geom_linerange(aes(ymin=p_value_lower, ymax=p_value_upper)) +
  labs(title="Simulation results for type I error with 5000 iterations using p value", x="sample size", y="type-I error")
ggsave(plot=plot, "figures/significance_p_value.jpg")

plot <- ggplot(res, aes(sample_size, power, color=method)) + 
  geom_line() + 
  geom_linerange(aes(ymin=lower, ymax=upper)) +
  labs(title="Simulation results for type I error with 5000 iterations using rsimsum", x="sample size", y="type-I error")
ggsave(plot=plot, "figures/significance_rsimsum.jpg")
