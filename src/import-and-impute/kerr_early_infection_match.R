## this is a quick script to highlight three important points:
## 1. there exist many people included in the analysis with symptom onset before the 7th of July
## 2. Most of the people that became ill before the 7th and then hospitalized were counted as non
## users
## 3. Hospitalised patient symptom onset data shows that the data is biased towards earlier
## infections *within* the study period
## written by Robin Mills
## Revised by Greg Tucker-Kellogg

library(here) #independent of Rstudio, sets here() to the git root
library(tidyverse)
library(glue)
library(janitor,quietly = TRUE)
source(here("src/data/global_params.R"), local = TRUE)
source(here("src/data/keep_results.R"), local = TRUE)

kerr_early_infection_tbl <- local({
    source(here("src/data/global_params.R"), local = TRUE)

                                        # let's start with what Kerr counts as infected participants
                                        # *and* non-participants. above 18 years old
    source(here("src/import-and-impute/kerr_infected.R"), local = TRUE)

                                        # Now get what Kerr counts as *all* IVM users kerr_users, regardless
                                        # of infection (data set 1)
    source(here("src/import-and-impute/kerr_users.R"), local = TRUE)

                                        # let's label the infected folks based on their use of IVM in the program
    kerr_infected$ivm_user <- kerr_infected$id %in% kerr_users$id


                                        # Now cross-check with Brazilian state data download this data from
                                        # google drive, too big for git
              
    source(here("src/import-and-impute", "INFLUD20.R"), local = TRUE)
    
    ## Allow for matches with city city-level cases, to cover non-residents who participated in the program
    ## Includes all DATASUS records from Itajai residents, reported by Itajai, or through Itajai hospitals in 2020.
    ## restricted only to records with hospitalisation or death records

    match_death_universe <- 
        citywide_infected |>
        filter(death) |> select(-hospitalised,-race,-type_2_diabetes)

    message(sprintf("There are %d records to consider for matching deaths between KC22 and DATASUS",nrow(match_death_universe)))
    results[['match_death_universe']] = nrow(match_death_universe)
    kerr_matched_deaths <-
        inner_join(match_death_universe,kerr_infected,
                   by = c("date_birth", "sex", "death"),
                   multiple = "all") |>
        filter(date_onset <= paper$end) |>
        arrange(desc(date_onset)) |>
        group_by(date_birth, sex, death) |>
        mutate(n_matches = n()) |>
        filter(row_number() == 1) |>
        ungroup() |>
        mutate(exposure = factor(ivm_user, labels = c("non-user", "user")),
               early = (date_onset < paper$start),
               early = factor(early, labels = c("false", "true"))) |>
        select(-ivm_user)

    message(sprintf("Matched %d death records between KC22 and DATASUS",nrow(kerr_matched_deaths)))
    results[['kerr_matched_deaths']] = nrow(kerr_matched_deaths)
    
    match_hosp_universe <- 
        citywide_infected |>
        filter(hospitalised,!death) |> select(-death,-race,-type_2_diabetes)

    message(sprintf("There are %d records to consider for matching hospitalised survivors between KC22 and DATASUS",nrow(match_hosp_universe)))
    results[['match_hosp_universe']] = nrow(match_hosp_universe)

    kerr_matched_hospitalisations <-
        inner_join(match_hosp_universe,kerr_infected,
                   by = c("date_birth", "sex", "hospitalised"),
                   multiple = "all") |>
        filter(date_onset <= paper$end) |>
        arrange(desc(date_onset)) |>
        group_by(date_birth, sex, hospitalised, death) |>
        group_by(date_birth, sex, death) |>
        mutate(n_matches = n()) |>
        filter(row_number() == 1) |>
        ungroup() |>
        mutate(exposure = factor(ivm_user, labels = c("non-user", "user")),
               early = (date_onset < paper$start),
               early = factor(early, labels = c("false", "true"))) |>
        select(-ivm_user)

    message(sprintf("Matched %d hospitalised survivor records between KC22 and DATASUS",nrow(kerr_matched_hospitalisations)))
    results[['kerr_matched_hospitalisations']] = nrow(kerr_matched_hospitalisations)
    
    ## This is the match between KC22 and the DATASUS data. We match on DOB, sex, hospitalisation, and
    ## death. In the case of multiple matches (could be two people, could be one person infected twice)
    ## we take the *latest* date of symptom onset to be conservative.

    results[['kerr_matched_cases']] = nrow(kerr_matched_hospitalisations) + nrow(kerr_matched_deaths)
    matched_hosp_and_death <- bind_rows(kerr_matched_deaths,kerr_matched_hospitalisations)
})

