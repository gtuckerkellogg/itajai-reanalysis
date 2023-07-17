#!/usr/bin/env/Rscript

library(here)
library(tidyverse)
library(gt)
source(here("src/functions","reporting.R"))

if (!exists("sim1")) load('results/sim1.RData')
if (!exists("sim2")) load('results/sim2.RData')
if (!exists("sim3")) load('results/sim3.RData')

ls()


sim1_infected_hosp <- map(sim1,filter,infected) |>
    map_df(report_hospitalisations) |>
    mutate(rate=true/(true+false)) |>
    group_by(exposure) |>
    summarise(model='i-ENR',outcome='hospitalisation',risk_ratio=median(estimate),
              upper=median(upper),
              lower=median(lower),
              p_value=median(fisher_exact),rate=median(rate))

sim1_infected_death <- map(sim1,filter,infected) |>
    map_df(report_deaths) |>
    mutate(rate=true/(true+false)) |>
    group_by(exposure) |>
    summarise(model='i-ENR',outcome='death',risk_ratio=median(estimate),
              upper=median(upper),
              lower=median(lower),
              p_value=median(fisher_exact),rate=median(rate))


sim3_infected_hosp <- map(sim3,filter,infected) |>
    map_df(report_hospitalisations,cutoff_date=as.Date("2020-09-30")) |>
    mutate(rate=true/(true+false)) |>
    group_by(exposure) |>
    summarise(model='i-KC22',,outcome='hospitalisation',risk_ratio=median(estimate),
              upper=median(upper),
              lower=median(lower),
              p_value=median(fisher_exact),rate=median(rate))

sim3_infected_death <- map(sim3,filter,infected) |>
    map_df(report_deaths,cutoff_date=as.Date("2020-09-30")) |>
    mutate(rate=true/(true+false)) |>
    group_by(exposure) |>
    summarise(model='i-KC22',outcome='death',risk_ratio=median(estimate),
              upper=median(upper),
              lower=median(lower),
              p_value=median(fisher_exact),rate=median(rate))


table_3 <- bind_rows(sim1_infected_hosp,
          sim1_infected_death,
          sim3_infected_hosp,          
          sim3_infected_death) |>
    mutate(CI=ifelse(is.na(upper),"",sprintf("(%3.2f -- %3.2f)",lower,upper))) |>
    select(-c(upper,lower)) |>
    relocate(CI,.after=risk_ratio)  |>
    gt(groupname_col='model',rowname_col='exposure') |>
    fmt_number(columns=c('risk_ratio')) |>
    fmt_percent(columns='rate') |>
    fmt_number(columns='p_value',n_sigfig = 2) |>
    sub_small_vals(columns=p_value,threshold=0.001) |>
    sub_missing(missing_text="")

gtsave(table_3,filename=here('results/table3.tex'))
    
