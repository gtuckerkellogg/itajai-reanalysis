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
    
    model_detailed <- citywide_infected |>  filter(date_notified >= paper$start,
                                                   date_onset >= paper$start - months(1),
                                                   date_onset <= paper$end
                                                   )

    model_imputed <- imputed_infections |> filter(date_notified >= paper$start,
                                                  date_onset >= paper$start -months(1),
                                                  date_onset <= paper$end
                                                  )

    create_infection_model(model_imputed,model_detailed)
})


##' simulate the Kerr et al study using notification dates recorded no earlier than July 7 and
##' death/hospitalisation sample from deaths and hospitalisiations
##' Uses a smoothed estimate of infection events from actual mapped deaths and hospitalisations
##' @title iKC22 model
##' @return function for cohort simulation.
##' @author Greg Tucker-Kellogg
sim3_model <- local({

    source(here("src/data/global_params.R"),local=TRUE)
    source(here("src/import-and-impute/simulated_enrolment.R"),local=TRUE)
    source(here("src/import-and-impute/kerr_early_infection_match.R"),local=TRUE)
    source(here("src/functions/event_simulation.R"),local=TRUE)
    source(here("src/functions/weibull_delay_generator.R"),local=TRUE)
    notification_delay <- delay_generator(dat_all,"date_onset","date_notified")
    
    
    model$start_date <- as.Date(paper$start - months(1))
    model$end_date <- paper$end
    model$study_period <-  interval(model$start_date,model$end_date)

    restrict_dates <- function(df) {
        filter(df,
               date_notified > paper$start,
               date_onset %within% model$study_period) |>
            arrange(date_onset)
        }

    ## this is to match Kerr
    model_detailed <- select(kerr_early_infection_tbl,
                             names(citywide_infected))

    late_events <- citywide_infected |>
        filter(ID_MN_RESI=="ITAJAI", ID_RG_RESI=="ITAJAI",
        ((date_hospitalised > model$end_date) | date_death > model$end_date),
        date_onset %within% paper$study_period)

    ## create a smoothed estimate of mapped deaths
    ## this needs to go further out
    mapped <- model_detailed |>
        tabyl(date_onset) |> as_tibble() |>
        mutate(days =as.numeric(date_onset - paper$start))

    model_days <-tibble(days=as.numeric(seq(model$start_date - months(1),model$end_date + months(1), by='days') - paper$start))
    
    mapped <- left_join(model_days,mapped) |> 
        mutate(n=replace_na(n,0),
               frac=n/sum(n),
               date_onset = paper$start+days)

    ## now predict using a loess smoother
    mapped$pred <- predict(loess(n ~ days,mapped,span=0.1,degree=0),model=TRUE);
    plot(mapped$date_onset,mapped$pred)
    mapped$.weights <- (mapped$pred - min(mapped$pred))/(max(mapped$pred) - min(mapped$pred))
    #plot(mapped$date_onset,mapped$.weights)
    overall_imputed <- imputed_infections |> restrict_dates()
    mapped_imputed <- mapped |> 
        slice_sample(n=nrow(imputed_infections),replace=TRUE,weight_by=.weights) |>
        mutate(date_notified=date_onset + notification_delay(nrow(imputed_infections))) |>
        select(date_notified,date_onset) |>
        restrict_dates()

    model_imputed <- mapped_imputed
    #model_imputed <- bind_rows(mapped_imputed,overall_imputed)

    model_detailed <- bind_rows(model_detailed,late_events)

    create_infection_model(model_imputed,model_detailed)
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
