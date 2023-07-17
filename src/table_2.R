#!/usr/bin/env Rscript
# Generate table 2 (the latex needs to be formatted, but it's there)

library(here)
library(tidyverse)
library(gt)
library(gtsummary)
library(readxl)

source(here("src/functions","reporting.R"))


source(here("src/import-and-impute/kerr_infected.R"))
source(here("src/import-and-impute/simulated_enrolment.R"))
source(here('src/functions/kerr_rr.R'),local=TRUE)

if (!exists("sim1")) load(here('results/sim1.RData'))
if (!exists("sim2")) load(here('results/sim2.RData'))
if (!exists("sim3")) load(here('results/sim3.RData'))

kerr_data <- 
    bind_rows(kerr_rr(),
          kerr_rr(exposures=c("non-user","user"),outcome='hospitalisation',infected_only=TRUE),
          kerr_rr(exposures=c("non-user","user"),outcome='death',infected_only=TRUE)
) |>
    mutate(model='KC22') |>
    mutate(kerr_risk_reduction=1-kerr_estimate,
           ci=ifelse(is.na(kerr_lower),NA,sprintf("[%3.2f--%3.2f]",kerr_lower,kerr_upper))) 


kerr_data |>
    filter(outcome=='death')  |> 
    summarise(total=sum(kerr_total),events=sum(kerr_true),rate=100*events/total)

tbl2_data <- bind_rows(
    report_infections(sim1) |> mutate(model='i-ENR'),
    report_infections(sim2) |> mutate(model='i-INF'),
    report_infections(sim3) |> mutate(model='i-KC22'),
    report_hospitalisations(sim1,infected_only=TRUE) |> mutate(model='i-ENR'),
    report_hospitalisations(sim2,infected_only=TRUE) |> mutate(model='i-INF'),
    report_hospitalisations(sim3,infected_only=TRUE) |> mutate(model='i-KC22'),
    report_deaths(sim1,infected_only=TRUE) |> mutate(model='i-ENR'),
    report_deaths(sim2,infected_only=TRUE) |> mutate(model='i-INF'),
    report_deaths(sim3,infected_only=TRUE) |> mutate(model='i-KC22'))


tbl2_data |>
    mutate(outcome=forcats::as_factor(outcome)) |> 
    select(exposure,outcome,model,matches('simulated')) |> 
    mutate(simulated_risk_reduction=1-simulated_estimate,
           ci=ifelse(is.na(simulated_lower),NA,sprintf("[%3.2f--%3.2f]",simulated_lower,simulated_upper))) ->
    tbl2_data_with_ci

ci_gt <- tbl2_data_with_ci |>
    select(exposure,outcome,model,ci) |>
    rbind(select(kerr_data,exposure,outcome,model,ci)) |>
    filter(exposure=='user') |>
    pivot_wider(names_from=model,values_from=ci) |>
    gt()

cat(as_latex(ci_gt))

tbl2_kerr <-
    kerr_data |> as_tibble() |> 
    mutate(outcome=forcats::as_factor(outcome))  |> 
    select(-ci,-subset) |>
    pivot_longer(cols=matches('kerr'),names_to='stat') |> 
    mutate(stat=stringr::str_remove(stat,"^kerr_"))  |>
    filter(!(stat %in% c('lower','upper'))) |>
    filter(exposure=='user' | stat == "rate",!(stat %in% c("true","false")))    

tbl2_sims <- tbl2_data_with_ci |>
    select(-ci) |>
    pivot_longer(cols=matches('simulated'),names_to='stat') |>
    mutate(stat=stringr::str_remove(stat,"^simulated_")) |>
    filter(exposure=='user' | stat == "rate",!(stat %in% c("true","false")))

main_gt <- 
    rbind(tbl2_sims,tbl2_kerr) |>
    filter(!(stat %in% c('total','lower','upper'))) |>
    mutate(stat=forcats::as_factor(stat)) |>
    pivot_wider(names_from=model,values_from=value)  |> 
    arrange(outcome,exposure,stat) |>
    mutate(stat=as.character(stat)) |>
    gt(groupname_col='outcome',rowname_col = 'stat') |>
    fmt_percent(columns=everything(),rows=matches('rate'),decimals=2) |>
    fmt_percent(columns=everything(),rows=matches('reduction'),decimals=0) |>
    fmt_number(columns=everything(),rows=matches('p_value'),decimals=3) |>
    fmt_number(columns=everything(),rows=matches('estimate'),decimals=2)


cat(map_chr(list(main_gt,ci_gt),as_latex),file=here("results/table2.tex"))

