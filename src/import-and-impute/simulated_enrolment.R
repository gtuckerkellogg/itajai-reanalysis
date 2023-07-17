library(here) #independent of Rstudio, sets here() to the git root
library(tidyverse)
library(readr)
library(readxl)
library(gt)
library(gtsummary)
library(stringi)

simulated_enrolment <- local({
    source(here("src","data","global_params.R"),local=TRUE)
## get the registration estimates from Brazil gov't and news

    registrations <-
        read_tsv(here('data','registrations-from-news.tsv'),show_col_types = FALSE) |>
        arrange(date) |>
        filter(!is.na(registrants)) |> 
        mutate(new_participants=registrants-lag(registrants),
               start_date=lag(date)+1,
               end_date=date,
               type='est') |>
        dplyr::select(date,registrants,start_date,end_date,new_participants) |>
        na.omit()
    
    registrations$date[nrow(registrations)] <- model$last_registration



## Now create simulated individual registrations assuming each window is a uniform distribution

simulated_participants <- nest(registrations,data=c(date,start_date,end_date,new_participants)) %>% 
    mutate(reg_date=lapply(data,function(df) {
        possible_dates <- seq(df$start_date,df$end_date,by='day')
        n <- df$new_participants
        sample(possible_dates,n,replace=TRUE)}
        )) %>%
    unnest(cols=c(data,reg_date)) %>%
    mutate(exposure='user') %>%
    select(exposure,reg_date) %>% sample_n(paper$users) # just subsample the number of participants +18yo

## we know how many non-particpants there were in the eventual data set, just put them all
## registering on day 0

simulated_nonparticipants <- tibble(exposure='non-user',
                                    reg_date=rep(paper$start,paper$non_users))


## now create a simulated cohort with unique IDs.

simulated_enrolment <- rbind(simulated_participants,simulated_nonparticipants) %>%
    mutate(id=stri_rand_strings(nrow(simulated_participants)+nrow(simulated_nonparticipants),16)) %>%
    select(id,exposure,reg_date)
simulated_enrolment
})

