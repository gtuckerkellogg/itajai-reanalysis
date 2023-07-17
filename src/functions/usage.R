## some functions for calculating and classifying ivermectin usage in simulated Kerr-et-al cohorts


library(here)
library(tidyverse)

##' return a modified cohort with an intended_usage column, derifed from sampling the actual tablet usage (up to
##' 80tabs) as in the supplemental data of KC22. non-users have 0 intende usage.
##' .. content for \details{} ..
##' @title intend_usage
##' @param cohort 
##' @return cohort (modified)
##' @author Greg Tucker-Kellogg
intended_usage <- local({
    source(here("src/import-and-impute","kerr_adult_participants.R"),local=TRUE)

    tablets <- table(kerr_adult_participants$total_tablets)
    n_tabs=as.integer(names(tablets))
    n_pax=as.integer(tablets)
    function(cohort) {
        cohort$intended_usage <- sample(n_tabs,size=nrow(cohort),replace=TRUE,prob=n_pax)
        cohort$intended_usage[cohort$exposure == 'non-user'] <- 0 
        cohort
    }
})

##' return a modified cohort with actual usage as the min(intended_usage,max_usage), where max_usage
##' is calculated by the time between enrolment and infection or study event_simulation.R nd date
##' @title actual_usage
##' @param cohort, including intended_usage
##' @param p_inf_stop
##' @return cohort (modified with additional column)
##' @author Greg Tucker-Kellogg
actual_usage <- local({
    source(here("src/data/global_params.R"),local=TRUE)
    treat_days <- rep(0,paper$duration)
    treat_days[c(seq(1,paper$duration,by=15),seq(2,paper$duration,by=15))] <- 1
    usage_days <- cumsum(treat_days)
    daily_usage <- 3
    function(cohort,p_inf_stop=model$p_inf_stop) {
        infected <- which(cohort$infected)
        hospitalised <- which(cohort$hospitalised)
        cohort$date_stop <- as.Date(NA)
        stop_on_onset <- infected[which(as.logical(rbinom(infected,1,p_inf_stop[cohort$intended_usage[infected]])))] 
        cohort$date_stop[stop_on_onset] <- cohort$date_onset[stop_on_onset]
        cohort$date_stop <- pmin(cohort$date_stop,cohort$date_hospitalised,cohort$date_death,model$end_date,na.rm=TRUE)
        duration <- pmax(1,as.numeric(interval(pmin(cohort$date_onset,cohort$reg_date,na.rm=TRUE),cohort$date_stop),'days'),na.rm=FALSE)
        max_usage <- usage_days[duration]*daily_usage
        cohort$actual_usage <- pmin(cohort$intended_usage,max_usage,na.rm=TRUE)
        cohort
    }
})
##' classify regularity by actual usage
##' @param cohort 
##' @return modified cohort
##' @author Greg Tucker-Kellogg
classify_regularity <- function(cohort) {
    mutate(cohort,
           regularity=ifelse(exposure=='non-user',"non-user",
                      ifelse(actual_usage>=30,"regular",
                      ifelse(actual_usage<=10,"irregular","other"))),
           regularity=fct_relevel(regularity,c('non-user','irregular','other','regular')))
    }
