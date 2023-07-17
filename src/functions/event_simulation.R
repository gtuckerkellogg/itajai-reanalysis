## This sets up a  common simulation framework

library(lubridate)
library(dplyr)
library(janitor)

##' Simulate (by sampling) the Kerr et al study using Itajai death data.
##' Generate a tibble with deaths sampled from Itajai matching the kerr count, and symtpom onset occurring the
##' study period. The number of observed deaths in the study period matches Kerr et al, but uncounted deaths
##' may occur after the study period.
##' @title death_simulation_model
##' @param model_detailed a tibble
##' @param cutoff_date 
##' @return a function to simulate deaths via a model
##' @author Greg Tucker-Kellogg (idea by Robin Mills)
death_simulation_model <- function(model_detailed,cutoff_date=paper$end) { 
    source(here("src/data/global_params.R"),local=TRUE)    
    death_table <- model_detailed |>
        filter(death,
               if_all(c(date_death,date_onset), ~ !is.na(.)))  |>
        mutate(id=NA,infected=TRUE) |> 
        select(id,infected,hospitalised,death,date_onset,date_hospitalised,date_death)        
    function() {
        death_table |> 
            slice_sample(n=2*paper$deaths_total,replace=TRUE) |>
            mutate(deaths_in_period=cumsum(date_death <= cutoff_date)) |>
            filter(deaths_in_period <= paper$deaths_total) |>
            select(-deaths_in_period)
    }
}



##' Simulate (by sampling) hospitalisations and deaths in the Kerr et al study by sampling from Itajai data.
##' Here the total number of hospitalisations during the study period matches Kerr et al, but uncounted hospitalisations
##' may also occur after the study period.
##' @title hosp_simulation_model
##' @param model_detailed a tibble
##' @param cutoff_date 
##' @return a function to simulate hospitalisations and deaths under a model
##' @author Greg Tucker-Kellogg (idea by Robin Mills)
hosp_simulation_model <- function(model_detailed,cutoff_date=paper$end) { 
    source(here("src/data/global_params.R"),local=TRUE)
    simulate_deaths <- death_simulation_model(model_detailed,cutoff_date)
    hospitalised_survivors <- model_detailed |>
        filter(hospitalised,!death,
               if_all(c(date_hospitalised,date_onset), ~ !is.na(.))) |>
        mutate(id=NA,infected=TRUE,hospitalised=TRUE,death=FALSE,date_death=NA) |> 
        select(id,infected,hospitalised,death,date_onset,date_hospitalised,date_death)        
    function() {
        repeat {
            death_table <- simulate_deaths()
            h_counted <- sum(death_table$date_hospitalised <= cutoff_date,na.rm=TRUE)
            to_count <- paper$hosp_total - h_counted
            to_sample <- 4*to_count
            if (to_count > 0) { break }
            }
        hosp_table <- hospitalised_survivors |> 
                             slice_sample(n=to_sample,replace=TRUE) |>
                             mutate(hosp_in_period=cumsum(date_hospitalised <=cutoff_date)) |>
                             filter(hosp_in_period <= to_count) |>
                             select(-hosp_in_period) |>
            bind_rows(death_table)
        hosp_table
    }
}
##' @title create_infection_model
##' @param imputed_infected a tibble
##' @param model_detailed a tibble
##' @param cutoff_date 
##' @return a function to simulate infections, hospitalisations, and deaths under a model
##' @author Greg Tucker-Kellogg
create_infection_model <-  function(model_imputed,model_detailed,cutoff_date=paper$end) { 
    source(here("src/data/global_params.R"),local=TRUE)    
    source(here("src/functions/usage.R"),local=TRUE)
    source(here("src/import-and-impute/simulated_enrolment.R"),local=TRUE)    
    w1 <- model_imputed %>% tabyl(date_onset)
    simulate_hospitalisations <- hosp_simulation_model(model_detailed,cutoff_date)
    function() {
        cohort <- mutate(simulated_enrolment,infected=FALSE,date_onset=as.Date(NA),
                         hospitalised=FALSE,date_hospitalised=as.Date(NA),
                         death=FALSE,date_death=as.Date(NA))

        ## the allocation function has a side effect of modifying cohort to prevent duplicate
        ## sampling. Takes a data frame shaped like the cohort,
        ## cohort members sampled accordingly
        
        allocate_date <- function(d,.n) {
            r <- sample(which((!cohort[['infected']]) & ((cohort$exposure == 'non-user' ) | (cohort$reg_date <= d))),.n)
            cohort$infected[r] <<- TRUE
            cohort$date_onset[r] <<- d
            cohort$id[r]
        }

        allocate <- function(the_tbl) {
            the_tbl %>%
                group_by(date_onset) %>%
                nest() %>% 
                mutate(data = map(data,
                                  ~ .x %>% mutate(id=allocate_date(date_onset,n())))) %>%
                unnest(cols=c(data)) %>%
                ungroup() %>%
                inner_join(select(cohort,id,exposure,reg_date),by='id') %>%
                select(colnames(cohort))
        }
        
        ## simulate hospitalisations and deaths and allocate to individuals in the cohort
        hospitalisations <- simulate_hospitalisations() %>% allocate()
        

        ## simulate non-hospitalised new infections and allocate
        new_infections <-
            tibble(id=NA,infected=TRUE,
                   date_onset=sample(w1$date_onset,
                                     size=(paper$covid_total-nrow(hospitalisations)),
                                     replace=TRUE,
                                     prob=w1$n),
                   hospitalised=FALSE,
                   death=FALSE,
                   date_hospitalised=as.Date(NA),
                   date_death=as.Date(NA)) |>
            allocate()
        
        infected <- rbind(new_infections,hospitalisations)

        rbind(infected,anti_join(cohort,infected,by='id')) |>
            intended_usage() |>
            actual_usage() |>
            classify_regularity()
    }
}





