#/usr/bin/env Rscript

## results showing immortal time bias in KC22 simulations even without the "regularity" analysis
## Only needs the 1000 simulations of each model 
library(here)
library(ggthemes)
library(patchwork)
library(gtsummary)
source(here("src/functions","reporting.R"))

if (!exists("sim1")) load('results/sim1.RData')


### propotion late

late_h <- sim1 %>%
    map_df(late_hospitalisation) %>% 
    mutate(proportion=(true/(true+false)),outcome='hospitalisation')

late_d  <- sim1 %>%
    map_df(late_death) %>%
    mutate(proportion=(true/(true+false)),outcome='death')

## get the proportion table

late_h %>% 
    select(exposure,false,true,proportion) %>%
    mutate(proportion=proportion*100) %>% 
    rename(in_time='false',late='true')  %>%
    tbl_summary(by='exposure') %>%
    add_p(include=proportion,test= proportion ~ 't.test')


late_d %>% 
    select(exposure,false,true,proportion) %>%
    mutate(proportion=proportion*100) %>% 
    rename(in_time='false',late='true')  %>%
    tbl_summary(by='exposure') %>%
    add_p(include=proportion,test= proportion ~ 't.test')


## get the risk ratio t test

late_d %>% filter(exposure=='user',is.finite(estimate)) %>% 
    pull(estimate) %>% t.test(mu=1)

late_h %>% filter(exposure=='user',is.finite(estimate)) %>% 
    pull(estimate) %>% t.test(mu=1)

late_h %>%
    filter(exposure=='user') %>%
    select(estimate) %>% 
    tbl_summary()



late_d %>%
    filter(exposure=='user') %>%
    select(estimate) %>% 
    tbl_summary()


proportion_late <- rbind(late_h,late_d) %>%
    mutate(outcome=fct_relevel(outcome,"hospitalisation")) %>% 
    ggplot(aes(x=outcome,y=proportion,fill=exposure)) + 
    scale_fill_colorblind() +
    geom_boxplot(alpha=0.5) +
    scale_y_continuous(labels = scales::percent,limits=c(0,0.4)) +
cowplot::theme_half_open() +
    theme(legend.position=c(0.05,.95),
          legend.title=element_blank()) + 
    ylab("Outcome events\nlost to attrition") + xlab("outcome")
proportion_late

rr_h <- sim1 %>%
    map_df(late_hospitalisation) %>%
    mutate(outcome='hospitalisation') %>%
    filter(exposure=='user') %>%
    select(outcome,estimate)

rr_d <- sim1 %>%
    map_df(late_death) %>%
    mutate(outcome='death') %>%
    filter(exposure=='user') %>%
    select(outcome,estimate)

late_risk <- rbind(rr_h,rr_d) %>%
    mutate(outcome=fct_relevel(outcome,"hospitalisation")) %>% 
    ggplot(aes(x=outcome,y=estimate,fill=outcome)) + 
    geom_hline(aes(yintercept=1),color='gray',linewidth=1) + 
    geom_boxplot(fill=colorblind_pal()(2)[2],alpha=0.5) +
    cowplot::theme_half_open() +
    scale_y_log10() +
    theme(legend.title=element_blank()) + 
    ylab("Relative risk\nof attition") + xlab("outcome")
late_risk


### first let's get a density plot 
date_onset_density <- sim1 %>%
    map_df(filter,infected) %>%
    ggplot(aes(x=date_onset,color=exposure,fill=exposure)) + geom_density(alpha=0.5,linewidth=2) +
    scale_fill_colorblind() +
    scale_color_colorblind() + 
    scale_y_continuous(labels=scales::percent_format(accuracy=1),
                       breaks=c(0,0.01,0.02)
                       ) + 
    facet_grid(exposure ~ .) + cowplot::theme_minimal_grid() +
    theme(legend.position ='none',
          strip.text.y=element_text(angle=90)) + 
    xlab("Onset date") + ylab("Density of\ninfection onsets") 
    
date_onset_density

# 3
example_death <- sim1[[3]] %>%
    filter(death) %>%
    group_by(exposure) %>%
    arrange(date_death,date_onset) %>%
    mutate(i=row_number(),r=i/n()) %>%
    select(i,exposure,date_onset,date_death) %>% 
    mutate(id=as.numeric(date_onset),
           hd=as.numeric(date_death)) %>%
    ggplot(aes(y=i,xmin=date_onset,xmax=date_death,color=exposure)) +
    geom_vline(aes(xintercept=paper$end),color='gray') + cowplot::theme_half_open() +
    geom_linerange() + 
    geom_point(aes(x=date_onset)) + 
    geom_point(aes(x=date_death)) +
    scale_x_date(date_breaks = "1 month", date_labels =  "%b")      + 
    xlab("Infection to death") +
    scale_color_colorblind() + 
    theme(axis.line.y = element_blank(),axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),axis.title.y=element_blank(),
          axis.text.x=element_text(angle=0, hjust=1),
          legend.position=c(-0.2,0.9),
          legend.title=element_blank())
example_death


time_bias_1 <- date_onset_density + example_death + proportion_late + late_risk + 
    plot_annotation(tag_levels="A")
time_bias_1

ggsave(here("results","attrition-bias.pdf"),
       plot=time_bias_1,
       width=8,
       height=6,
       units="in")


