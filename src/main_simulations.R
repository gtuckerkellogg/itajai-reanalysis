#!/usr/bin/env Rscript

### This is effectively a generative resampling of simulated infections 
## set seed to the date of Ivermectin approval for human use in the US ;-)
## https://www.thepharmaletter.com/article/usa-clears-ivermectin-for-human-use

library(here)
library(tidyverse)
library(tictoc)
set.seed(as.numeric(as.Date("1996-02-12")))
source(here("src/simulation_models.R"),local=TRUE)

iterations <- 10

## This is the iENR model

tic(sprintf("sim1 (%d times)",iterations),quiet=FALSE)
sim1 <- n_rounds_sim(iterations,sim1_model)
toc()
tic(sprintf("sim1 (iENR) saving"))
save(sim1,file=here('results/sim1.RData'))
toc()
rm(sim1)

## This runs the iINF model

tic(sprintf("sim2 (%d times)",iterations))
sim2 <- n_rounds_sim(iterations,sim2_model)
toc()
tic(sprintf("sim2 (iINF) saving"))
save(sim2,file=here('results/sim2.RData'))
rm(sim2)
toc()


## This runs the iKC22 model

tic(sprintf("sim3 (%d times)",iterations))
sim3 <- n_rounds_sim(iterations,sim3_model)
toc()
tic(sprintf("sim3 (iKC22) saving"))
save(sim3,file=here('results/sim3.RData'))
rm(sim3)
toc()

