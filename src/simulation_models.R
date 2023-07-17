## Three simulation models with different setups for infection dates
## All using Itajai-only data
## Author: Greg Tucker-Kellogg

library(here)
library(tidyverse)
library(lubridate)
library(janitor)
library(progressr)
handlers(global=TRUE)


# This must be sourced first because citywide_infected and itajai_cases are used to create locally specific
# versions in each of the models.
source(here("src/import-and-impute","INFLUD20.R"),local=TRUE)
source(here("src/import-and-impute","itajai_cases.R"),local=TRUE)        

##' simulate the Kerr et al study using symptom onset dates recorded *within* the study period
##' @title sim1_model
##' @return simulation of a cohort
##' @author Greg Tucker-Kellogg
sim1_model <- local({
    source(here("src/data/global_params.R"),local=TRUE)
    source(here("src/import-and-impute/simulated_enrolment.R"),local=TRUE)
    source(here("src/functions/event_simulation.R"),local=TRUE)

    model_detailed <- citywide_infected |> filter(date_onset %within% paper$study_period)
    model_imputed <- imputed_infections |> filter(date_onset %within% paper$study_period) 
   create_infection_model(model_imputed,model_detailed)
})

##' simulate the Kerr et al study using notification dates recorded within the study period. Symtpom onset
##' dates may precede the study period.
##' @title sim2_model
##' @return simulation of a cohort
##' @author Greg Tucker-Kellogg
sim2_model <- local({
    source(here("src/data/global_params.R"),local=TRUE)
    source(here("src/import-and-impute/simulated_enrolment.R"),local=TRUE)
    source(here("src/functions/event_simulation.R"),local=TRUE)
    
    model_detailed <- citywide_infected |>  filter(date_notified %within% paper$study_period,
                                                   date_onset >= paper$start - months(1))

    model_imputed <- imputed_infections |> filter(date_notified %within% paper$study_period,
                                                  date_onset >= paper$start -months(1))

    create_infection_model(model_imputed,model_detailed)
})

##' simulate the Kerr et al study using notification dates recorded no earlier than July 7 and all
##' observations stopped on September 30.
##' Symptom onset dates may precede the study period.
##' @title sim3_model
##' @return simulation of a cohort
##' @author Greg Tucker-Kellogg
sim3_model <- local({
    source(here("src/data/global_params.R"),local=TRUE)
    source(here("src/import-and-impute/simulated_enrolment.R"),local=TRUE)
    source(here("src/import-and-impute/kerr_early_infection_match.R"),local=TRUE)
    source(here("src/functions/event_simulation.R"),local=TRUE)
    
    model$start_date <- as.Date(paper$start - months(1))
    model$end_date <- ymd("2020-09-30")
    model$study_period <- interval(model$start_date,model$end_date)
    ## this is to match Kerr
    model_detailed <- kerr_early_infection_tbl |> filter(date_onset %within% model$study_period)
    model_imputed <- imputed_infections |> filter(date_notified %within% paper$study_period,
                                                  date_onset %within% model$study_period)

                                        #function() simulate_infections()
    create_infection_model(model_imputed,model_detailed,model$end_date)
})
##' siumulate the Kerr et al study using smoothed weighted probability of infeciton from citywide records
##' weighted by the empirical rate of inclusion in KC22
##'
##' .. content for \details{} ..
##' @title iKC22_model
##' @return simulation of a cohort
##' @author Greg Tucker-Kellogg
iKC22_model <- local({

    source(here("src/data/global_params.R"),local=TRUE)
    source(here("src/import-and-impute/simulated_enrolment.R"),local=TRUE)
    source(here("src/import-and-impute/kerr_early_infection_match.R"),local=TRUE)
    source(here("src/functions/event_simulation.R"),local=TRUE)
    
    model$start_date <- as.Date(paper$start - months(1))
    model$end_date <- ymd("2020-09-30")
    model$study_period <- interval(model$start_date,model$end_date)
    model_detailed <- kerr_early_infection_tbl |> 
        filter(date_notified %within% paper$study_period,
               date_onset >= paper$start - months(1)) |> 
        arrange(date_onset)


    ## create a smoothed estimate of mapped hospitalisations
    ## this needs to go further out

    mapped <-
        kerr_early_infection_tbl |>
        dplyr::select(date_onset,date_notified) |>
        tabyl(date_onset) |> as_tibble() |>
        mutate(days=as.numeric(date_onset-paper$start))

    days <- tibble(days=seq(min(mapped$days),max(mapped$days),1))
    mapped <- left_join(days,mapped) |>
        mutate(n=replace_na(n,0),
               frac=n/sum(n),
               date_onset=paper$start+days)

    ## here's the loess smoother

    mapped.lo <- loess(n ~ days,mapped,span=0.3,degree=2)

    ## now predict
    mapped_pred <- predict(mapped.lo,data.frame(days=days),se=TRUE) |>
        as_tibble() |> 
        bind_cols(days) |>
        mutate(fit=ifelse(fit<0,0,fit),
               date_onset=paper$start+days) |>
        select(date_onset,fit)

    ggplot(mapped_pred,aes(date_onset,fit)) + geom_line()

    impute_weights  <- citywide_infected |>
        filter(date_notified %within% paper$study_period) |>
#        filter(hospitalised) |> 
        tabyl(date_onset) |>
        left_join(mapped_pred) |>
        as_tibble() |> 
        filter(date_onset %within% model$study_period) |>
        replace_na(list(fit=0)) |>
        mutate(weight=fit/n) |>
        select(date_onset,weight) |>
        left_join(tabyl(citywide_infected,date_onset)) |>
        mutate(fit=weight*n) 

    ggplot(impute_weights,aes(date_onset,fit)) + geom_line()

    model_imputed <- tibble(date_onset=sample(impute_weights$date_onset,nrow(imputed_infections),
                                              replace=TRUE,prob=impute_weights$fit))
    create_infection_model(model_imputed,model_detailed,model$end_date)
})


##' Generate multiple simulations of the Kerr et al study using specified models
##' @title n_round_sim
##' @param n: the number of rounds of simulation
##' @param simulation_model: a function for generating a single simulation under a given model framework
##' @return a list of simulated cohorts
n_rounds_sim  <- local({
##' @author Greg Tucker-Kellogg
      function(n,simulation_model) {
        xs <- 1:n
        p <- progressor(along=xs)
        map(xs,function(x) {
            p(sprintf("x=%g", x))
            simulation_model()
        })
    }
})

rm(dat_all,citywide_infected,imputed_infections,itajai_cases)
