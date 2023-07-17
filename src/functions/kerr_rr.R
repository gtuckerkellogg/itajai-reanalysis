#!/usr/bin/env Rscrpt

library(here)
library(tidyverse)


kerr_rr <- local({

    kerr_summary <- readr::read_tsv(here('data/kerr/table-summary.tsv')) |>
        tidyr::replace_na(list(subset="all",outcome="")) |>
        select(-`exposure level`) |> 
        mutate(exposure=forcats::as_factor(exposure),
               infected_only=ifelse(subset=="infected",TRUE,FALSE),
               psm=ifelse(subset=="infected (PSM)",TRUE,FALSE))
    
    function(exposures=c('non-user','user'),outcome='infection',infected_only=FALSE) {

        .outcome <- outcome
        .infected_only <- infected_only

        totals <-
            kerr_summary |>
            filter(exposure %in% exposures,outcome=="",infected_only==.infected_only,!psm)  |>
            select(exposure,total=n)

        event_table <- 
            kerr_summary |> 
            filter(exposure %in% exposures,outcome==.outcome,infected_only==.infected_only,!psm) |>
            select(exposure,subset,true=n) |>
            inner_join(totals,by=join_by(exposure)) |>
            mutate(false=total-true,rate=true/total) |> 
            select(false,true) |> 
            as.matrix()

        RR <- epitools::riskratio.boot(event_table)
        cohort <-
            cbind(totals,RR$data[1:2,1:2],as_tibble(RR$p.value)[,'fisher.exact']) |>
            as_tibble() |>
            mutate(rate=true/total) |>
            rename(p_value='fisher.exact')

        cbind(cohort,as_tibble(RR$measure)) |>
            rename_with(~paste0("kerr_",.x),everything()) |>
            mutate(outcome=.outcome,subset=ifelse(.infected_only,"infected","all")) |>
            rename(exposure='kerr_exposure') |>
            relocate(exposure,outcome,subset)
    }
})
