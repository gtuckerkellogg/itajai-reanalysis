## Import and clean itajai_cases from Brazil national data,
## and impute infection and onset dates

library(here) 
library(tidyverse)
library(lubridate)

## actual case information from SUS
## onset dates imputed

if (!exists("dat_all")) source(here("src/import-and-impute","INFLUD20.R"),local=TRUE)

itajai_cases <- local({
    source(here("src/data/global_params.R"),local=TRUE)
    source(here("src/functions/weibull_delay_generator.R"),local=TRUE)
    itajai <- here("data",'itajai_june_to_december.csv')

    colspec <- spec_delim(itajai,delim=";")
    colspec$cols[['date']] <- col_date(format="%d/%m/%Y")
    itajai_cases <- read_delim(itajai,
                               col_types=colspec) %>%
        select(date,new_cases,new_deaths)
    rm(itajai,colspec)

    ## The raw data includes negative cases. Shift the negative cases to the day before
    for (i in which(itajai_cases$new_cases < 0)) {
        itajai_cases$new_cases[i-1] <-     itajai_cases$new_cases[i-1] -  itajai_cases$new_cases[i]
        itajai_cases$new_cases[i] <- 0 
    }

    ## There is a huge backlog of cases on 31 August. 
    ## The backlog is so large it'll create an artifact  let's get the average of the surrounding
    ## two months and distribute accordingly


    target <- itajai_cases %>%
        filter(date > "2020-07-30",date < "2020-10-01", date != "2020-08-31") %>% 
        pull(new_cases) %>%
        mean() %>%
        round()

    ## 47 cases

    excess <- max(itajai_cases$new_cases) - target
    idx <- last(order(itajai_cases$new_cases))
    days_to_distribute <- 60

    itajai_cases$new_cases[idx] <- target
    itajai_cases$new_cases[(idx - days_to_distribute):(idx-1)] <- itajai_cases$new_cases[(idx-days_to_distribute):(idx-1)] + table(sample(1:days_to_distribute,excess,replace=TRUE))

    select(itajai_cases,date,new_cases,new_deaths)
})

imputed_infections <- local({
    source(here("src/data/global_params.R"),local=TRUE)    
    source(here("src/functions/weibull_delay_generator.R"),local=TRUE)
    notification_delay <-delay_generator(dat_all,"date_onset","date_notified")
    tibble(date_notified=do.call("c",map2(itajai_cases$date,itajai_cases$new_cases, function(d,n) rep(d,n))),
           date_onset=do.call("c",map2(itajai_cases$date,itajai_cases$new_cases,function(d,n)  { as.Date(d) - notification_delay(n)})))
})




