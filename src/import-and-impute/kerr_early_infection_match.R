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
library(janitor)
kerr_early_infection_tbl <- local({
    source(here("src/data/global_params.R"),local=TRUE)

                                        # let's start with what Kerr counts as infected participants *and*
                                        # non-participants. above 18 years old
    source(here('src/import-and-impute/kerr_infected.R'),local=TRUE)

                                        # Now get what Kerr counts as *all* IVM users kerr_users, regardless
                                        # of infection (data set 1)
    source(here('src/import-and-impute/kerr_users.R'),local=TRUE)

                                        # let's label the infected folks based on their use of IVM in the program
    kerr_infected$ivm_user <- kerr_infected$id %in% kerr_users$id


                                        # Now cross-check with Brazilian state data download this data from
                                        # google drive, too big for git
    if (!exists("citywide_infected")) source(here("src/import-and-impute","INFLUD20.R"),local=TRUE)


    ## This is the match between KC22 and the DATASUS data. We match on DOB, sex, hospitalisation, and
    ## death. In the case of multiple matches (could be two people, could be one person infected twice)
    ## we take the *latest* date of symptom onset to be conservative. 

    kerr_infected |>
        inner_join(citywide_infected,by=c('date_birth','sex','hospitalised','death'), multiple="all")  |>
        filter(date_onset <= paper$end) |> 
        arrange(desc(date_onset)) |>
                                        # If multiple matches, conservatively choose the latest possible
                                        # matched date of onset
        group_by(date_birth,sex,hospitalised,death) |>
        mutate(n_matches=n()) |>
        filter(row_number()==1) |>
        ungroup() |>
        mutate(exposure=factor(ivm_user,labels=c("non-user","user")),
               early=(date_onset < paper$start),
               early=factor(early,labels=c('false','true'))) |>
        select(-ivm_user)

})

