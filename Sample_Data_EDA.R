# Set working directory ---------------------------------------------------
rm(list=ls())
library(tidyverse)
library(rms)
setwd("F:/yifei/UCL stats for clinical trial/Simulation with Rafael/Simulations")


# load the df
df <- read_csv("data/sample_data.csv")

# Plot histogram of all outcomes
ggplot(data=df, aes(x=outcome)) + geom_histogram() + scale_x_continuous(limits=c(0,40), breaks = seq(0, 40, by=5))
histo_by_trt <- ggplot(data=df, aes(x=outcome)) + 
  geom_histogram() + 
  scale_x_continuous(limits=c(0,40), breaks = seq(0, 40, by=5)) + 
  facet_wrap(~treatment) +
  labs(title="Histogram of simulated outcome by treatment arm, size 20000")
histo_by_trt

# Frequency polygon of outcomes by treatment
freq_poly <- ggplot(data=df, aes(x=outcome, color=factor(treatment))) + 
  geom_freqpoly(binwidth=1) + 
  scale_x_continuous(limits=c(0,40), breaks = seq(0, 40, by=5)) + 
  theme(legend.position=c(.9, .75)) +
  labs(title="Freqeuncy polygon of simulated outcome, size 20000")
freq_poly

# Plot histograms of logged values
df <- df %>% 
  mutate(log_outcome = log(outcome))
log_histo_by_trt <- ggplot(data=df, aes(x=log_outcome)) + 
  geom_histogram() + 
  facet_wrap(~treatment) +
  labs(title="Histogram of logged simulated outcome by treatment arm, size 20000")
log_histo_by_trt

# Look at the summary statistics for outcome by treatment arm
df %>% 
  group_by(treatment) %>% 
  summarise(prop_zero=mean(outcome==0.01), 
            mean=mean(outcome), 
            median=median(outcome),
            sd=sd(outcome),
            q_25=quantile(outcome, 0.25),
            q_75=quantile(outcome, 0.75))

# Test code for using linear regression
# treatment coefficients: mean difference
model <- lm(outcome ~ treatment, data=df)
summ <- summary(model)
print(summ)
#print(paste0("expected: ", mean_intervention - mean_control))
#p_value <- summ$coefficients[2,4]
#print(p_value)

# Test code for using linear regression on log-transformed outcome
# e^intercept = control median
# e^(intercept+treatment) = intervention median
# treatment coefficient: log ratio of median
df$log_outcome <- log(df$outcome)
model <- lm(log_outcome ~ treatment, data=df)
summ <- summary(model)
print(summ)
#p_value <- summ$coefficients[2,4]

# Test code for using GLM (gaussian with log link)
# e^intercept = control mean
# e^(intercept+treatment) = intervention mean
# treatment coefficient: log ratio of mean
model <- glm(outcome ~ treatment, data=df, family=gaussian(link="log"))
summ <- summary(model)
print(summ)
#p_value <- summ$coefficients[2,4]
#print(p_value)

# Test code for using GLM (Gamma with log link)
# e^intercept = control mean
# e^(intercept+treatment) = intervention mean
# treatment coefficient: log ratio of mean
model <- glm(outcome ~ treatment, data=df, family=Gamma(link="log"))
summ <- summary(model)
print(summ)
#p_value <- summ$coefficients[2,4]
#print(p_value)

# Test code for using ordinal logistic regression
# ??Interpretation of coefficients??
# log odds??
dd <- datadist(df)
options(datadist = "dd")
model <- orm(outcome ~ treatment, data=df)
print(model)
