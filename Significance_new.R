# Set working directory ---------------------------------------------------
rm(list=ls())
library(tidyverse)
library(nleqslv)
library(glue)
library(rms)
library(rsimsum)
setwd("F:/yifei/UCL stats for clinical trial/Simulation with Rafael/Simulations")

# get a list of all files in functions folder
pathnames <- list.files(pattern="[.]R$", path="function", full.names=TRUE)
# source all of the files in functions folder
sapply(pathnames, FUN=source)

result <- run_simulation(
  control_median=6,
  control_variance=16,
  control_zero_prop=0.14,
  intervention_median=6,
  intervention_variance=16,
  intervention_zero_prop=0.14,
  num_iter=10000,
  sample_size_start=100,
  sample_size_end=160,
  sample_size_step=5,
  function_paths=pathnames
)


write_csv(result, glue("F:/yifei/UCL stats for clinical trial/Simulation with Rafael/Simulations/result/significance_new.csv"))
