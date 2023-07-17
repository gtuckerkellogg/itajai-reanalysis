#!/usr/bin/env Rscript
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
library(readr)
library(readxl)
library(gt)
library(viridis)
library(colorspace)
library(ggthemes)
library(colorblindr)
source(here("src/data/global_params.R"), local=TRUE)
library(gtsummary)
library(patchwork)
source(here("src/import-and-impute", "INFLUD20.R"), local=TRUE)
source(here("src", "import-and-impute", "kerr_early_infection_match.R"), local=TRUE)
source(here("src/import-and-impute/simulated_enrolment.R"))


table1a <-  kerr_early_infection_tbl |>
     select(early,exposure) |>
     mutate(early=fct_recode(early,'Covid onset on or after 7 July 2020'='false',
                             'Covid onset before 7 July 2020'='true')) |>
    tbl_summary(by=exposure,label=list(early ~ "All Brazilian Health Ministry-mapped individuals")) |>
    modify_header(label ~ "") |>
    bold_labels() |>
    add_p(~ "fisher.test") |> as_gt()

# table1a


gtsave(table1a,here('results','early_infections.tex'))

table1b <- kerr_early_infection_tbl |>
    select(early,exposure,death) |>
    mutate(death=as_factor(death)) |>
    mutate(death=fct_recode(death,'Alive'='FALSE',
                             'Dead'='TRUE')) |>
    mutate(early=fct_recode(early,'Covid onset on or after 7 July 2020'='false',
                            'Covid onset before 7 July 2020'='true')) |>
    tbl_strata(strata=death,
               .tbl_fun= ~ .x |>
                   tbl_summary(by=exposure,
                               label=list(early ~ "")) |>
                   modify_header(label ~ "",all_stat_cols() ~ "**{level}**") |>
                   bold_labels() |>                   
                   add_p(~ "fisher.test"),
               .combine_with='tbl_stack',
               ) |>
    as_gt()

 table1b

gtsave(table1b,here('results','early_infections_deaths.tex'))


x_start <- as.Date("2020-04-01")
x_end <- as.Date("2020-12-31")
common_xlim <- xlim(x_start,x_end)

.colors <- colorblind_pal()(2)
study_rect <- annotate("rect",xmin=paper$start,xmax=paper$end,ymin=-Inf,ymax=Inf,fill=.colors[2],alpha=0.2)
bias_rect_1 <- annotate("rect",xmin=as.Date("2020-04-01"),xmax=paper$start,ymin=-Inf,ymax=Inf,
                      fill=.colors[1],alpha=0.2)
bias_rect_2 <- annotate("rect",xmin=paper$start,xmax=as.Date("2020-09-30"),ymin=-Inf,ymax=Inf,
                      fill=.colors[2],alpha=0.2)


matched_p <- kerr_early_infection_tbl |>
    filter(!(hospitalised & is.na(date_hospitalised))) |>
    ggplot(aes(x=date_onset)) +
    study_rect +
    geom_histogram(binwidth=2,alpha=0.5,fill=.colors[1]) +
    annotate("text",x=as.Date("2020-09-19"),y=15,label="KC22 study period") + 
    xlab("Symptom onset date") + ylab("hospitalised adults") +
    common_xlim  + 
    ggtitle("Hospitalised Itajaí Covid-19 cases (matched to KC22)") +
    cowplot::theme_half_open() 
#matched_p <- cowplot::ggdraw(edit_colors(matched_p,tritan))
 matched_p


citywide_infected |> 
    filter(hospitalised,date_birth<=paper$start - years(18)) |>
    filter(date_onset >= as.Date("2020-01-01"),date_onset <= as.Date("2020-12-31")) |>
    filter(date_onset >= paper$start,date_onset <= paper$end)  |>
    summarise(deaths=sum(death),hospitalised=sum(hospitalised),death_rate_of_hosp=deaths/hospitalised)



    

citywide_p <- citywide_infected |>
    filter(hospitalised,date_birth<=paper$start - years(18))  |>
    ggplot(aes(x=date_onset)) +
    study_rect + 
    geom_histogram(binwidth=2,alpha=0.5,fill=.colors[1]) +
    annotate('text',x=as.Date("2020-09-19"),y=60,label="KC22 study period") + 
    common_xlim +
    xlab("Symptom onset date") + 
    ylab("hospitalised adults") +
    ggtitle("Hospitalised Itajaí Covid-19 cases (DATASUS)") + 
    cowplot::theme_half_open()

### now let's get a probability of enrolment

enrolment <- simulated_enrolment |>
    arrange(reg_date)  |>
    group_by(reg_date,exposure) |>
    summarise(daily_enrolment=n()) |>
    ungroup() |>
    group_by(exposure) |> 
    mutate(enrolment=cumsum(daily_enrolment)) |> 
    pivot_wider(names_from=exposure,values_from=enrolment,values_fill=0) |> 
    select(-daily_enrolment) |> 
    mutate(`non-user`=paper$non_users) 

enrolment_p <- (head(enrolment,1) |> mutate(reg_date=as.Date("2020-04-01"))) |> 
    rbind(enrolment) |>
    rbind(tail(enrolment,1) |> mutate(reg_date=paper$end)) |> 
    janitor::clean_names() |>
    rowwise() |> 
    mutate(p_nonuser=non_user/(non_user+user)) |> 
    ggplot(aes(x=reg_date,y=p_nonuser)) + 
    bias_rect_1 +
    bias_rect_2 + 
    geom_line() +
    annotate("text",y=0.1,x=(paper$start - (paper$start - as.Date("2020-04-01"))/2),label="Prior infections") + 
    annotate("text",x=(paper$start + (as.Date("2020-09-30") - paper$start)/2),y=0.1,label="Enrolment bias") + 
    common_xlim +
    xlab("Enrolment date") +
    ylab("p(non-user)") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1),limits=c(0,1)) + 
    ggtitle("Probability of non-user group allocation") + 
    cowplot::theme_half_open()

citywide_p + matched_p + enrolment_p + 
    plot_layout(nrow=3)  +
    plot_annotation(tag_levels = 'A')

ggsave(here("results","early_infections.pdf"),
       width=8,
       height=7,
       units="in")


