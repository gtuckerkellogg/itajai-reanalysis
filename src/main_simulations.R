#!/usr/bin/env Rscript

### This is effectively a generative resampling of simulated infections 
## set seed to the date of Ivermectin approval for human use in the US ;-)
## https://www.thepharmaletter.com/article/usa-clears-ivermectin-for-human-use

library(here)
library(tidyverse)
library(tictoc)
library(glue)
set.seed(as.numeric(as.Date("1996-02-12")))
source(here("src/simulation_models.R"),local=TRUE)

gc <- function() {
    invisible(base::gc())
}

iterations <- 1000

models <- c("sim1", "sim2", "sim3")

for (model_name in models) {
    print(tic(glue("{model_name}: ({iterations} times)"),quiet=FALSE))
    sim_model <- eval(sym(glue("{model_name}_model")))
    cohorts <- n_rounds_sim(iterations,sim_model)
    toc()
    tic(glue("saving {model_name}"))
    write_rds(cohorts,file=here("results",glue("{model_name}.rds")))
    toc()
}

