#!/usr/bin/env Rscript
## this is a quick script to highlight three important points:
## 1. there exist many people included in the analysis with symptom onset before the 7th of July
## 2. Most of the people that became ill before the 7th and then hospitalized were counted as non
## users
## 3. Hospitalised patient symptom onset data shows that the data is biased towards earlier
## infections *within* the study period
## written by Robin Mills
## Revised by Greg Tucker-Kellogg

suppressPackageStartupMessages({
library(here) # independent of Rstudio, sets here() to the git root
library(tidyverse)
library(janitor)
library(readr)
library(readxl)
library(gt)
library(viridis)
library(colorspace)
library(ggthemes)
library(scales)
library(colorblindr)
source(here("src/data/global_params.R"), local = TRUE)
source(here("src/data/keep_results.R"), local = TRUE)
library(gtsummary)
library(patchwork)
options(dplyr.summarise.inform = FALSE)
})
source(here("src/import-and-impute", "INFLUD20.R"), local = TRUE)
source(here("src", "import-and-impute", "kerr_early_infection_match.R"), local = TRUE)
source(here("src/import-and-impute/simulated_enrolment.R"))

write_tsv(kerr_early_infection_tbl,file=here("results/kerr_sus_match.tsv"))
  
message("Comparing in 2nd half 2020")

onset_in_2h <- citywide_infected |>
    filter(ID_MN_RESI=="ITAJAI", ID_RG_RESI=="ITAJAI") |>
    filter(date_birth <= paper$start - years(18)) |>
    filter(date_onset >= as.Date("2020-07-01"), date_onset <= as.Date("2020-12-31"))

results[['deaths_after_onset_in_2H2020']] <- sum(onset_in_2h$death)
results[['hosp_after_onset_in_2H2020']] <- sum(onset_in_2h$hospitalised)

deaths_in_2h <- citywide_infected |>
    filter(ID_MN_RESI=="ITAJAI", ID_RG_RESI=="ITAJAI") |>
    filter(date_birth <= paper$start - years(18)) |>
    filter(death,date_death >= as.Date("2020-07-01"), date_death <= as.Date("2020-12-31"))

results[['deaths_2H_2020']] <- sum(deaths_in_2h$death)
results[['hosp_2H_2020']] <- sum(deaths_in_2h$hospitalised)

message("Comparing in study period")

onset_in_period <- citywide_infected |>
    filter(ID_MN_RESI=="ITAJAI", ID_RG_RESI=="ITAJAI") |>
    filter(date_birth <= paper$start - years(18)) |>
    filter(date_onset %within% paper$study_period)

onset_in_period |>
    summarise(deaths = sum(death), hospitalised = sum(hospitalised), death_rate_of_hosp = deaths / hospitalised)

hospitalised_in_period <- citywide_infected |>
    filter(ID_MN_RESI=="ITAJAI", ID_RG_RESI=="ITAJAI") |>
    filter(date_birth <= paper$start - years(18)) |>
    filter(hospitalised,date_hospitalised %within% paper$study_period)
deaths_in_period <- citywide_infected |>
    filter(ID_MN_RESI=="ITAJAI", ID_RG_RESI=="ITAJAI") |>
    filter(date_birth <= paper$start - years(18)) |>
    filter(death,date_death %within% paper$study_period)


results[['deaths_after_onset_in_study_period']] <- sum(onset_in_period$death)
results[['hosp_after_onset_in_study_period']] <- sum(onset_in_period$hospitalised)
results[['deaths_in_study_period']] <- nrow(deaths_in_period)
results[['hosp_in_study_period']] <- nrow(hospitalised_in_period)

missed_deaths <-
    anti_join(onset_in_period,kerr_early_infection_tbl,by=c('date_birth', 'sex', 'death')) |>
    filter(death)

results[['missed_deaths']] <- nrow(missed_deaths)

table1a <- kerr_early_infection_tbl |>
    select(early, exposure) |>
    mutate(early = fct_recode(early,
                              "Covid onset on or after 7 July 2020" = "false",
                              "Covid onset before 7 July 2020" = "true"
                              )) |>
    tbl_summary(
        by = exposure,
        label = list(early ~ "All Brazilian Health Ministry-mapped individuals")
    ) |>
    modify_header(label ~ "") |>
    bold_labels() |>
    add_p(~"fisher.test") |>
    as_gt()

#table1a


gtsave(table1a, here("results", "early_infections.tex"))

table1b <- kerr_early_infection_tbl |>
    select(early, exposure, death) |>
    mutate(death = as_factor(death)) |>
    mutate(death = fct_recode(death,
                              "Alive" = "FALSE",
                              "Dead" = "TRUE"
                              )) |>
    mutate(early = fct_recode(early,
                              "Covid onset on or after 7 July 2020" = "false",
                              "Covid onset before 7 July 2020" = "true"
                              )) |>
    tbl_strata(
        strata = death,
        .tbl_fun = ~ .x |>
            tbl_summary(
                by = exposure,
                label = list(early ~ "")
            ) |>
            modify_header(label ~ "", all_stat_cols() ~ "**{level}**") |>
            bold_labels() |>
            add_p(~"fisher.test"),
        .combine_with = "tbl_stack",
        ) |>
    as_gt()

#table1b

gtsave(table1b, here("results", "early_infections_deaths.tex"))


x_start <- as.Date("2020-04-01")
x_end <- as.Date("2020-12-31")
common_xlim <- xlim(x_start, x_end)

.colors <- colorblind_pal()(3)
study_rect <- annotate("rect", xmin = paper$start, xmax = paper$end, ymin = -Inf, ymax = Inf, fill = .colors[2], alpha = 0.2)
bias_rect_1 <- annotate("rect",
                        xmin = as.Date("2020-04-01"), xmax = paper$start, ymin = -Inf, ymax = Inf,
                        fill = .colors[1], alpha = 0.2
                        )
bias_rect_2 <- annotate("rect",
                        xmin = paper$start, xmax = as.Date("2020-09-30"), ymin = -Inf, ymax = Inf,
                        fill = .colors[2], alpha = 0.2
                        )


ann_study <- tibble(date = as.Date("2020-09-19"), y = 6.5, lab ="KC22 study period",event=factor("death",levels=c("onset","hospitalised","death")))
#ann_study

event <- factor(c("onset","hospitalised","death"))
annots <- tibble(date = as.Date("2020-12-05"), y = c(8,6.6,3),event=event,lab=c("Symptom\nonset","Hospitalisation","Death"))

matched_p <- kerr_early_infection_tbl |>
    pivot_longer(cols=c(date_onset,date_hospitalised,date_death),
                 names_to='event',values_to = 'date') |> select(event,date,death) |>
    mutate(event=as_factor(str_replace(event,"date_","")))  |> 
    ggplot(aes(x = date)) +
    study_rect +
    ylab("count") + xlab("Event date") + 
    geom_histogram(binwidth = 2,na.rm=TRUE) + scale_fill_hc() +
    facet_wrap(fct_rev(event) ~ .,nrow=3,scales='free_y',strip.position='right',dir='h') +
    scale_y_continuous(breaks= pretty_breaks()) +
    cowplot::theme_half_open() + 
    ggtitle("Event dates of matched Itajaí COVID-19 hospitalisations and deaths") +
    common_xlim + 
    theme(strip.background = element_blank(),strip.placement = 'inside',
          strip.text = element_blank()) + 
    geom_text(data=ann_study,aes(y=y,label=lab)) +
    geom_text(data=annots,aes(y=y,label=lab),hjust=0.0,size=6) 

#matched_p

fig_s1 <- missed_deaths |>
    pivot_longer(cols=c(date_onset,date_hospitalised,date_death),
                 names_to='event',values_to = 'date') |> select(event,date,death) |>
    mutate(event=as_factor(str_replace(event,"date_","")))  |> 
    ggplot(aes(x = date)) +
    study_rect +
    annotate("text", x = as.Date("2020-09-19"), y = 12, label = "KC22 study period") +    
    ylab("count") + xlab("Event date") +
    geom_histogram(binwidth = 2,na.rm=TRUE) + scale_fill_hc() +
    facet_wrap(event ~ .,nrow=3,scales='free_y',strip.position='left') +
    scale_y_continuous(breaks= pretty_breaks()) +
    cowplot::theme_half_open() +
    common_xlim +
    theme(strip.background = element_blank(),
          strip.placement = 'inside',
          strip.text.y.left = element_text(angle = 0,hjust=0,vjust=1)) 


ggsave(here("results", "figure_s1.pdf"),
       fig_s1,
       width = 7,
       height = 5,
       units = "in"
       )

citywide_p <- citywide_infected |>
    filter(hospitalised | death, date_birth <= paper$start - years(18)) |>
    ggplot(aes(x = date_onset)) +
    study_rect +
    annotate("text", x = as.Date("2020-09-19"), y = 36, label = "KC22 study period") +
    xlab("Symptom onset date") +
    ylab("count") +
    geom_histogram(binwidth = 2, na.rm=TRUE,alpha = 0.6, fill = .colors[1]) +
    common_xlim +
    ggtitle("Symptom onset date for Itajaí cases leading to hospitalisation or death (DATASUS)") +
    cowplot::theme_half_open()

                                        #citywide_p

### now let's get a probability of enrolment

enrolment <- simulated_enrolment |>
    arrange(reg_date) |>
    group_by(reg_date, exposure) |>
    summarise(daily_enrolment = n()) |>
    ungroup() |>
    group_by(exposure) |>
    mutate(enrolment = cumsum(daily_enrolment)) |>
    pivot_wider(names_from = exposure, values_from = enrolment, values_fill = 0) |>
    select(-daily_enrolment) |>
    mutate(`non-user` = paper$non_users)

enrolment_p <- (head(enrolment, 1) |> mutate(reg_date = as.Date("2020-04-01"))) |>
    rbind(enrolment) |>
    rbind(tail(enrolment, 1) |> mutate(reg_date = paper$end)) |>
    janitor::clean_names() |>
    rowwise() |>
    mutate(p_nonuser = non_user / (non_user + user)) |>
    ggplot(aes(x = reg_date, y = p_nonuser)) +
    bias_rect_1 +
    bias_rect_2 +
    geom_line() +
    annotate("text", y = 0.1, x = (paper$start - (paper$start - as.Date("2020-04-01")) / 2), label = "Prior infections") +
    annotate("text", x = (paper$start + (as.Date("2020-09-30") - paper$start) / 2), y = 0.1, label = "Enrolment bias") +
    common_xlim +
    xlab("Enrolment date") +
    ylab("p(non-user)") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(0, 1)) +
    ggtitle("Probability of non-user group allocation") +
    cowplot::theme_half_open()

figure_1 <- matched_p + citywide_p + enrolment_p +
    plot_layout(nrow = 3,heights=c(3.5,1,1),guides='collect') + 
    plot_annotation(tag_levels = "A")
#figure_1

ggsave(here("results", "early_infections.pdf"),
       figure_1,
       width = 15,
       height = 12,
       units = "in"
       )
