#!/usr/bin/env Rscript
### This is effectively a generative resampling of simulated infections 

library(here)
library(tidyverse)
source(here("src/data/global_params.R"))
source(here("src","functions","reporting.R"))
source(here("src","functions","usage.R"))

library(ggthemes)
library(patchwork)
library(gtsummary)
library(washu)
library(glue)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
   model <- "sim1"
} else if (length(args)==1) {
    tryCatch({model <- match.arg(args[1],c("sim1","sim2","sim3"))},
             error=function(e) e)
}

stopifnot(model %in% c("sim1","sim2","sim3"))

message(glue("******************   model is {model}   ******************"))

if (!exists(model)) cohorts <- read_rds(glue('results/{model}.rds'))

## make sure the regularity is correct.

sim_data <- cohorts |>
    map(mutate,actual_usage=intended_usage) |>
    map(actual_usage) |>
    map(classify_regularity)

reg_death_report <- sim_data %>%
    map_df(report_deaths_cohort,exposures=c("non-user","regular")) %>%
    mutate(proportion=true/(true+false)) |>
    as_tibble()

reg_hosp_report <- sim_data %>%
    map_df(report_hospitalisations_cohort,exposures=c("regular","non-user")) %>%
    mutate(proportion=true/(true+false)) |>
    as_tibble()

versus_nonuser <- rbind(reg_hosp_report,reg_death_report) %>%
    filter(exposure=='regular') %>%
    mutate(comparison="versus non-user",
           outcome=fct_relevel(outcome,"hospitalised"),
           outcome=fct_recode(outcome,hospitalisation="hospitalised"))

reg_vs_irreg_death <- sim_data %>% 
    map_df(report_deaths_cohort,exposures=c("regular","irregular")) %>%
    mutate(proportion=true/(true+false))

reg_vs_irreg_hosp <- sim_data %>% 
    map_df(report_hospitalisations_cohort,exposures=c("regular","irregular")) %>%
    mutate(proportion=true/(true+false))

versus_irregular <- rbind(reg_vs_irreg_death,reg_vs_irreg_hosp) %>%
    filter(exposure=='regular') %>%
    mutate(comparison="versus irregular users",
           outcome=fct_relevel(outcome,"hospitalised"),
           outcome=fct_recode(outcome,hospitalisation="hospitalised"))


## irreg_v_nonuser_death <- sim_data %>% 
##     map_df(report_deaths_cohort,exposures=c("irregular","non-user")) %>%
##     mutate(proportion=true/(true+false))

## irreg_v_nonuser_hosp <- sim_data %>% 
##     map_df(report_hospitalisations_cohort,exposures=c("irregular","non-user")) %>%
##     mutate(proportion=true/(true+false))

## irreg_vs_nonuser <- rbind(irreg_v_nonuser_death,irreg_v_nonuser_hosp) |>
##     filter(exposure=='irregular') %>%
##     mutate(comparison="versus irregular users",
##            outcome=fct_relevel(outcome,"hospitalised"),
##            outcome=fct_recode(outcome,hospitalisation="hospitalised"))


reg_color <- colorblind_pal()(2)[2]

reguser_protection <- rbind(versus_nonuser,versus_irregular) |>
    mutate(outcome=fct_recode(outcome,hospital="hospitalisation")) %>% 
    ggplot(aes(x=outcome,y=estimate,fill=exposure)) + 
    scale_fill_colorblind() + 
    geom_hline(aes(yintercept=1),color='gray',linewidth=1)  + 
    facet_grid(~ comparison) + 
    geom_boxplot(fill=reg_color) +
    scale_y_log10() +        
    cowplot::theme_half_open()  + 
    theme(legend.position='none',
          legend.title=element_blank()) +
    ylab("Relative risk of outcome") + xlab("outcome") 

hd_by_group <- map_df(sim_data,filter,hospitalised) %>% 
    mutate(regularity=fct_relevel(regularity,'non-user')) %>%
    ggplot(aes(x=date_hospitalised,y=regularity,fill=regularity)) +
    geom_vline(aes(xintercept=paper$end),color='gray') +
    ylab("") + xlab("date of hospitalisation") + 
    geom_boxplot(fill=colorblind_pal()(4)[c(1,3,4,2)]) + 
    cowplot::theme_half_open() 


dd_by_group <- map_df(sim_data,filter,death) %>% 
    mutate(regularity=fct_relevel(regularity,'non-user')) %>%
    ggplot(aes(x=date_death,y=regularity,fill=regularity)) +
    geom_vline(aes(xintercept=paper$end),color='gray') +
    ylab("") + xlab("date of death") + 
    geom_boxplot(fill=colorblind_pal()(4)[c(1,3,4,2)]) + 
    cowplot::theme_half_open()

reguser_protection_p <- ((hd_by_group / dd_by_group) | reguser_protection) +
    plot_annotation("Immortal time bias protection of infected 'regular' users",
                    theme = theme(plot.title = element_text(size = rel(1.5))),
                    tag_levels = "A")

ggsave(glue(here("results","reguser_protection"),"-{model}.pdf"),
       plot = reguser_protection_p,
       width=9,
       height=5,
       units="in")

