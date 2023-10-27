#!/usr/bin/env Rscript
# Generate table 2 (the latex needs to be formatted, but it's there)

library(here)
library(tidyverse)
library(gt)
library(gtsummary)
library(readxl)
library(glue)

source(here("src/functions","reporting.R"))


source(here("src/import-and-impute/kerr_infected.R"))
source(here("src/import-and-impute/simulated_enrolment.R"))
source(here('src/functions/kerr_rr.R'),local=TRUE)


table2_latex <- function(gt_obj) {
    as_latex(gt_obj) |>
        as.character() |>
        str_replace_all(c(
            "[$]" = "",
            "longtable" = "tabulary",
            "0.000" = "<0.001"
        )) 
    }


infected_only_p <- FALSE


model_names <- tibble(model=c('sim1','sim2','sim3'),
                      name=c('i-ENR', 'i-INF','i-KC22'))

tbl2_data <- map_df(1:nrow(model_names),
                    function(i) {
                        message(glue("running {model_names$model[i]}"))
                        cohorts <- read_rds(here("results", sprintf("%s.rds", model_names$model[i])))
                        bind_rows(report_infections(cohorts),
                                  report_hospitalisations(cohorts,infected_only=infected_only_p),
                                  report_deaths(cohorts,infected_only=infected_only_p)) |>
                            mutate(model=model_names$name[i])
                    }) |>
    mutate(outcome=as_factor(outcome),
           outcome=fct_relevel(outcome,"infection","hospitalisation"),
           ) |> arrange(outcome) 

sim_data <- tbl2_data |> 
    select(outcome,model,exposure,matches('simulated')) |>
    rename_with(~ str_remove(. , "simulated_"),matches("simulated_")) |>
    select(-true,-false)

kerr_data <- bind_rows(
    kerr_rr(exposures=c("non-user","user"),outcome='infection'),
    kerr_rr(exposures=c("non-user","user"),outcome='hospitalisation',infected_only=infected_only_p),
    kerr_rr(exposures=c("non-user","user"),outcome='death',infected_only=infected_only_p)) |>
    mutate(model='KC22') |>
    rename_with(~ str_remove(. , "kerr_"),matches("kerr_")) |>
    select(names(sim_data))
        
neg_data <- bind_rows(
    kerr_neg(exposures=c("non-user","user"),outcome='infection'),
    kerr_neg(exposures=c("non-user","user"),outcome='hospitalisation',infected_only=infected_only_p),
    kerr_neg(exposures=c("non-user","user"),outcome='death',infected_only=infected_only_p)) |>
    mutate(model='i-NEG',lower=NA,upper=NA,estimate=NA,p_value=NA) |>
    rename_with(~ str_remove(. , "kerr_"),matches("kerr_"))  |>
    select(names(sim_data))
        


ci_df <- bind_rows(sim_data,kerr_data) |>
    select(outcome,model,exposure,lower,upper) |>
    filter(exposure=='user') |>
    mutate(ci  = sprintf("[%3.2f--%3.2f]",lower,upper),stat='rate') |>
    select(-lower,-upper) |> 
    pivot_wider(names_from=model,values_from=ci)    

main_gt <- bind_rows(neg_data,sim_data,kerr_data) |>
    mutate(outcome=fct_relevel(as_factor(outcome),"infection","hospitalisation")) |>
    select(-lower,-upper) |>
    mutate(risk_reduction=1-estimate) |> 
    pivot_longer(cols=c(rate,estimate,p_value,risk_reduction),names_to='stat') |>
    filter(exposure=='user' | stat =='rate') |>
    mutate(stat=relevel(as_factor(stat),ref='rate')) |>
    arrange(stat) |> 
    pivot_wider(names_from=model,values_from=value) |>
    mutate(exposure=ifelse(stat=='rate',exposure,"")) |>
    gt(groupname_col='outcome',rowname_col = 'stat')  |>
    fmt_percent(columns=everything(),rows= stat == 'rate', decimals=2) |>
    fmt_percent(columns=everything(),rows= stat == 'risk_reduction', decimals=0) |>
    fmt_number(rows = stat %in% c('p_value', 'estimate'), decimals = 3) |>
    sub_small_vals(rows = stat == 'p_value', , threshold = 0.001)

cat(table2_latex(main_gt),as_latex(gt(ci_df)),file=here("results/table2.tex"))

