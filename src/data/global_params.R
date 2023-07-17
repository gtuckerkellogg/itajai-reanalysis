library(lubridate)



model <- list(
    ## this is the max delay (in days) that is included when generating Weibull distributions
    delay_limit=30,
    ## chance of stopping on infection
#    p_inf_stop=c(rep(0.5,29),rep(0.1,80-29)),
    p_inf_stop=c(rep(0.3,29),rep(0.05,80-29)),
    ## chance of stopping on hospitalisation
    p_hosp_stop=rep(1,80),
    seed=as.numeric(as.Date("1996-02-12")),
    last_registration=lubridate::ymd("2020-09-30")
)



paper = list()

paper$start <- as.Date("2020-07-07")
paper$end <- as.Date("2020-12-02")
model$end_date <- paper$end
paper$study_period <- interval(paper$start,paper$end)
paper$duration <- paper$study_period / days(1)

paper$users = 113844 
paper$non_users = 45716 
paper$regular_users = 8325
paper$irregular_users = 33971
paper$total_participants = paper$non_users + paper$users
paper$cases <- paper$total_participants
# covid statistics
paper$covid_users = 4197
paper$covid_non_users = 3034
paper$covid_irregular_users = 1542
paper$covid_regular_users = 283
paper$covid_total = 3034 + 4197 

# hospitalization statistics 
paper$hosp_users = 86
paper$hosp_non_users = 99
paper$hosp_irregular_users = 38
paper$hosp_regular_users = 0
paper$hosp_total = 185

# death statistics
paper$deaths_users = 62
paper$deaths_non_users = 79 
paper$deaths_irregular_users = 29
paper$deaths_regular_users = 2
paper$deaths_total = paper$deaths_users + paper$deaths_non_users


## city stats
paper$itajai_adults=161545
