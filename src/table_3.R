#!/usr/bin/env Rscript

library(here)
library(tidyverse)
library(gt)
library(glue)
source(here("src/functions","reporting.R"))



protection <- map_df(c("sim1", "sim2", "sim3"),
                     function(kerr_model) {
                         print(glue("reading {kerr_model}"))
                         cohorts <- read_rds(here(glue("results/{kerr_model}.rds")))
           bind_rows(
               mutate(late_hosp_protection(cohorts),model = kerr_model,outcome='hospitalisation'),
               mutate(late_death_protection(cohorts),model = kerr_model,outcome='death'))
                     })

table_3 <- protection |>
    relocate(model,outcome,exposure) |>
    mutate(model= as.character(model_names[model])) |>
    gt(groupname_col='model',rowname_col='exposure')  |>
    fmt_number(rows = everything(), columns = matches("estimate|lower|upper")) |>
    fmt_number(rows = everything(), columns = c("p","outcome_rate"), decimals = 3) |>
    fmt_percent(rows = everything(), columns = "outcome_rate", decimals = 0) |>
    sub_small_vals(rows = everything(), columns = "p", threshold = 0.001) |>
    cols_merge_range(
        col_begin = lower,
        col_end = upper)  |>
    cols_label(lower="95% CI",
               estimate = "Risk Ratio",
               outcome_rate = "rate of lost events")

gtsave(table_3,filename=here('results/table3.tex'))
