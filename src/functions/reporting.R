## simulate IPD infections in a cohort.
## smooth to get something continuous, simulate infections as a Poisson
## simulate delay times as Weibull

library(here)
library(tidyverse)
library(lubridate)
library(washu)
source(here("src/data/global_params.R"),local=TRUE)
source(here('src/functions/kerr_rr.R'),local=TRUE)
library(progressr)
handlers(global = TRUE)

clean_rr <- function(cohort,rr,outcome) {
    outcome <- rlang::ensym(outcome)
    dn <- list(unique(dplyr::pull(cohort, exposure)), unique(dplyr::pull(cohort,outcome)))
    names(dn) <- c("exposure",outcome)
    dplyr::bind_cols(tibble::tibble(exposure = sort(dn$exposure), outcome= rlang::as_string(outcome)),
                     tibble::as_tibble(rr$data)[-dim(rr$data)[1], 
                                                -dim(rr$data)[2]],
                     tibble::as_tibble(rr$measure), tibble::as_tibble(rr$p.value),
                     tibble::tibble(correction = rr$correction, method = attr(rr,"method"))) %>%
        janitor::clean_names()
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param cohort 
##' @param cutoff_date 
##' @return 
##' @author Greg Tucker-Kellogg
truncate_cohort <- function(cohort,cutoff_date=paper$end) {
    mutate(cohort,
           infected=(infected & (date_onset <= cutoff_date)), 
           hospitalised=(hospitalised & (date_hospitalised <= cutoff_date)),
           death=(death & (date_death <= cutoff_date)))
}

##' compute a corrected risk ratio, confidence interval, and p value from 
##' two contingency tables, one experimental and one simulated
##' Each individual risk ratio is calculated by bootstrap estimate, then the ratio of boostraps
##' is used to calculate a corrected RR with confidence interval and p value.
##' @title corrected estimate
##' @param M_e a 2x2 matrix of experimental observations. 
##' @param M_s a 2x2 matrix of experimental observations. 
##' @param conf.level default 95% 
##' @return a matrix of corrected estimates
##' @author Greg Tucker-Kellogg
corrected_estimate <- function(M_e,M_s,conf.level=0.95) {
    alpha <- 1 - conf.level
    rr_boot <- function(C_n,C_e,I_n,I_e,replicates = 10000) {
        C_tot <- C_n + C_e
        C_p <- C_e/C_tot
        I_tot <- I_n + I_e
        I_p <- I_e/I_tot
        r_i <- rbinom(replicates, I_tot, I_p)/I_tot
        r_c <- rbinom(replicates, C_tot, C_p)/C_tot
        rrboot <- r_i/r_c
        rrboot
    }
    rrboot_e <- rr_boot(M_e[1,1],M_e[1,2],M_e[2,1],M_e[2,2])
    rrboot_s <- rr_boot(M_s[1,1],M_s[1,2],M_s[2,1],M_s[2,2])
    rrboot_c <- rrboot_e/rrboot_s
    rr <- median(rrboot_c,na.rm=TRUE)
    p.val <- sum(rrboot_c >= 1)/length(rrboot_c)
    CI <- quantile(rrboot_c, c(alpha/2, 1 - alpha/2), na.rm = T)
    matrix(c(1,1,1,NA,rr,CI,p.val),byrow=TRUE,nrow=2,
           dimnames=list(NULL,paste("corrected",c("estimate","lower","upper","p"),sep="_")))
}


##' return the SD of the RR on a log scale provided a matrix of counts
##' @title 
##' @param m
##' @return a number, on a natural log scale, of the SD
##' @author Greg Tucker-Kellogg
LRR_SD <- function(m) {
    m <- as.matrix(m)
    as.numeric(sqrt(m[2,1]/(m[2,2]*(m[2,2]+m[2,1])) +
         m[1,1]/(m[1,2]*(m[1,2]+m[1,1]))))
}


report_infections <- function(cohorts,truncate=TRUE,cutoff_date=paper$end,exposures=c("non-user","user"))  {
    kerr <- kerr_rr(exposures=exposures)
    M_e <- select(kerr,kerr_false,kerr_true) |> as.matrix()    
    p <- progressor(along = 1:length(cohorts))
    report_cohort <- function(cohort,truncate=TRUE,cutoff_date=paper$end,exposures=c("non-user","user")) {
        p(sprintf("infections: (%s)",paste(exposures,collapse=" vs ")))
        if (!all(exposures %in% c("non-user","user"))) {
            cohort <- remap_regularity_to_exposure(cohort,exposures)
            }
        if(truncate) { 
            cohort <-truncate_cohort(cohort,cutoff_date)
        }
        rr <- epitools::riskratio.boot(cohort$exposure,cohort$infected,replicates=10000)
        M_s <- janitor::tabyl(cohort,exposure,infected) |> select("FALSE","TRUE") |> as.matrix()
        infections <- clean_rr(cohort,rr,infected) |>
            cbind(corrected_estimate(M_e,M_s))
        infections
    }
    

    infection_tbl <- map_df(cohorts,report_cohort,truncate=truncate,cutoff_date=cutoff_date,exposures=exposures) |>
        dplyr::mutate(rate=true/(true+false)) |>
        group_by(exposure) |>
        reframe(simulated_estimate=median(estimate),
                simulated_lower=median(lower),
                simulated_upper=median(upper),
                simulated_true=round(median(true)),
                simulated_false=round(median(false)),
                simulated_rate=median(rate),
                simulated_p_value=median(fisher_exact),
                corrected_estimate=median(corrected_estimate,na.rm=TRUE),
                corrected_lower=median(corrected_lower,na.rm=TRUE),
                corrected_upper=median(corrected_upper,na.rm=TRUE),
                corrected_p_value=median(corrected_p,na.rm=TRUE),
                outcome='infection') |>
        inner_join(kerr,by=join_by(exposure, outcome))
    infection_tbl$comparison <- paste(infection_tbl$exposure,collapse=' vs ')
    infection_tbl
}

report_hospitalisations <- function(cohorts,truncate=TRUE,cutoff_date=paper$end,infected_only=FALSE,exposures=c("user","non-user")) {
    kerr <- kerr_rr(exposures=exposures,outcome='hospitalisation',infected_only=infected_only)
    M_e <- select(kerr,kerr_false,kerr_true) |> as.matrix()    
    p <- progressor(along = 1:length(cohorts))
    report_cohort <- function(cohort) { #},truncate=TRUE,cutoff_date=paper$end,infected_only=infected_only,exposures=c("user","non-user")) {
        p(sprintf("hospitalisations: (%s)",paste(exposures,collapse=" vs ")))
        if (infected_only) {
            cohort <- filter(cohort,infected)
        }
        if (!all(exposures %in% c("non-user","user"))) {
            cohort <- remap_regularity_to_exposure(cohort,exposures)
            }
       if(truncate) { 
            cohort <-truncate_cohort(cohort,cutoff_date)
        }
        rr <- with(cohort,epitools::riskratio.boot(exposure,hospitalised,replicates=10000))
        M_s <- janitor::tabyl(cohort,exposure,hospitalised) |> select("FALSE","TRUE") |> as.matrix()        
        hospitalisations <- clean_rr(cohort,rr,hospitalised) |>
            cbind(corrected_estimate(M_e,M_s))
        hospitalisations
    }

    h_tbl <- map_df(cohorts,report_cohort) |> #,truncate=truncate,cutoff_date=cutoff_date,infected_only=infected_only,exposures=exposures) |>
        dplyr::mutate(rate=true/(true+false)) |> 
        group_by(exposure) |>
        summarise(simulated_estimate=median(estimate),
                  simulated_upper=median(upper),
                  simulated_lower=median(lower),
                  simulated_true=round(median(true)),
                  simulated_false=round(median(false)),
                  simulated_rate=median(rate),                  
                  simulated_p_value=median(fisher_exact),
                  corrected_estimate=median(corrected_estimate,na.rm=TRUE),
                  corrected_lower=median(corrected_lower,na.rm=TRUE),
                  corrected_upper=median(corrected_upper,na.rm=TRUE),
                  corrected_p_value=median(corrected_p,na.rm=TRUE),
                  outcome='hospitalisation') |> 
       inner_join(kerr,by=join_by(exposure, outcome))
    h_tbl$comparison <- paste(h_tbl$exposure,collapse=' vs ')
    h_tbl
}

report_deaths <- function(cohorts,truncate=TRUE,cutoff_date=paper$end,infected_only=FALSE,exposures=c("user","non-user")) {
    kerr <- kerr_rr(exposures=exposures,outcome='death',infected_only=infected_only)
    M_e <- select(kerr,kerr_false,kerr_true) |> as.matrix()
    p <- progressor(along = 1:length(cohorts))
    report_cohort <- function(cohort) { # truncate=TRUE,cutoff_date=paper$end,infected_only=TRUE,exposures=c("user","non-user")) {
        p(sprintf("deaths: (%s)",paste(exposures,collapse=" vs ")))
        if (infected_only) {
            cohort <- filter(cohort,infected)
        }
        if (!all(exposures %in% c("non-user","user"))) {
            cohort <- remap_regularity_to_exposure(cohort,exposures)
            }
       if(truncate) { 
            cohort <-truncate_cohort(cohort,cutoff_date)
        }
        rr <- with(cohort,epitools::riskratio.boot(exposure,death,replicates=10000))
        M_s <- janitor::tabyl(cohort,exposure,death) |> select("FALSE","TRUE") |> as.matrix()                
        deaths <- clean_rr(cohort,rr,death) |>
            cbind(corrected_estimate(M_e,M_s))
        deaths
    }

    d_tbl <- map_df(cohorts,report_cohort) |> #truncate=truncate,cutoff_date=cutoff_date,infected_only=infected_only,exposures=exposures) |>
        mutate(rate=true/(true+false)) |> 
        group_by(exposure) |>
        summarise(simulated_estimate=median(estimate),
                  simulated_upper=median(upper),
                  simulated_lower=median(lower),
                  simulated_true=round(median(true)),
                  simulated_false=round(median(false)),
                  simulated_rate=median(rate),                  
                  simulated_p_value=median(fisher_exact),
                  corrected_estimate=median(corrected_estimate,na.rm=TRUE),
                  corrected_lower=median(corrected_lower,na.rm=TRUE),
                  corrected_upper=median(corrected_upper,na.rm=TRUE),
                  corrected_p_value=median(corrected_p,na.rm=TRUE),
                  outcome='death') |> 
       inner_join(kerr,by=join_by(exposure, outcome))    
    d_tbl$comparison <- paste(d_tbl$exposure,collapse=' vs ')
    d_tbl

}


late_hospitalisation <- function(cohort,cutoff_date=paper$end) { 
    cohort %>%
        filter(hospitalised) %>%
        mutate(late_hosp=factor(as.character(date_hospitalised > cutoff_date),levels=c("FALSE","TRUE"))) -> 
        cohort
    or <- epitools::oddsratio.wald(cohort$exposure,cohort$late_hosp)
    clean_rr(cohort,or,"late_hosp")
}

late_death <- function(cohort) {
    cohort %>%
        filter(death) %>%
        mutate(late_death=factor(as.character(date_death > paper$end),levels=c("FALSE","TRUE"))) -> 
        cohort
    or <- epitools::oddsratio.wald(cohort$exposure,cohort$late_death)
    clean_rr(cohort,or,'late_death')
}


median_infection_dates <- function(cohort) {
    cohort %>% group_by(exposure) %>%
        summarise(median_idate=median(date_onset,na.rm=TRUE))
}

##' For reporting purposes, filter and remap regularity to exposure
##' 
##' .. content for \details{} ..
##' @title 
##' @param cohort 
##' @param exposures 
##' @param infected_only 
##' @return 
##' @author Greg Tucker-Kellogg
remap_regularity_to_exposure <- function(cohort,exposures=c("regular","non-user")) {
    filter(cohort,regularity %in% exposures) |> 
        select(-exposure) |> 
        mutate(exposure=droplevels(regularity))
    }

report_deaths_regularity <- function(cohort,truncate=TRUE,infected_only=TRUE,exposures=c("regular","non-user")) {
    if (truncate) {
        cohort <- mutate(cohort,death=(death & (date_death <= paper$end))) 
    }
    if (infected_only) { cohort <- filter(cohort,infected) }
    cohort <- filter(cohort,regularity %in% exposures) |>
        mutate(regularity=fct_relevel(droplevels(regularity),exposures))
    rr <- with(cohort,epitools::riskratio.boot(regularity,death,replicates=10000))
    clean_rr(cohort,rr,'death') |> mutate(exposure=exposures)
}

report_hospitalisations_regularity <- function(cohort,truncate=TRUE,infected_only=TRUE,exposures=c("regular","non-user")) {
    if (truncate) {
        cohort <- mutate(cohort,hospitalised=(hospitalised & (date_hospitalised <= paper$end))) 
    }
    if (infected_only) { cohort <- filter(cohort,infected) }    
    if (infected_only) { cohort <- filter(cohort,infected) }
    cohort <- filter(cohort,regularity %in% exposures) |>
        mutate(regularity=fct_relevel(droplevels(regularity),exposures))
    rr <- with(cohort,epitools::riskratio.boot(regularity,hospitalised,replicates=10000))
    clean_rr(cohort,rr,'hospitalised') %>% mutate(exposure=sort(exposures))
}

