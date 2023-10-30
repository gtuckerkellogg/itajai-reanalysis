#!/usr/bin/env Rscript

library(here)
library(tidyverse)
library(gt)
library(gtsummary)
source(here("src", "data/global_params.R"))
source(here("src/functions", "reporting.R"))
source(here("src/functions", "usage.R"))
library(progressr)
library(glue)
handlers(global = TRUE)

gc <- function() {
    invisible(base::gc())
    invisible(base::gc())
    invisible(base::gc())
}

 reassign_regularity <- function(cohorts, probabilistic = TRUE) {
    if (probabilistic) {
        P_INF_STOP <- model$p_inf_stop
    } else {
        P_INF_STOP <- model$p_hosp_stop
    }
    map(cohorts, actual_usage, p_inf_stop = P_INF_STOP) |>
        map(classify_regularity)
}

##' This generates a single data frame for all the relevant outcomes and exposures. It's meant to be turned
##' into a table for publication
##' @title overview_outcomes
##' @param cohorts a list of cohorts
##' @param cutoff_date the date after which evens should not be considered
##' @return a data frame of outcomes and risk ratios
##' @author Greg Tucker-Kellogg
overview_outcomes <- function(cohorts, cutoff_date = model$end, infected_only = FALSE) {
    gc()
    rbind(
        report_infections(cohorts, cutoff_date = cutoff_date, exposures = c("non-user", "user")),
        report_infections(cohorts, cutoff_date = cutoff_date, exposures = c("non-user", "irregular")),
        report_infections(cohorts, cutoff_date = cutoff_date, exposures = c("non-user", "regular")),
        report_infections(cohorts, cutoff_date = cutoff_date, exposures = c("irregular", "regular")),
        report_hospitalisations(cohorts, cutoff_date = cutoff_date, infected_only = infected_only, exposures = c("non-user", "user")),
        report_hospitalisations(cohorts, cutoff_date = cutoff_date, infected_only = infected_only, exposures = c("non-user", "irregular")),
        report_hospitalisations(cohorts, cutoff_date = cutoff_date, infected_only = infected_only, exposures = c("non-user", "regular")),
        report_hospitalisations(cohorts, cutoff_date = cutoff_date, infected_only = infected_only, exposures = c("irregular", "regular")),
        report_deaths(cohorts, cutoff_date = cutoff_date, infected_only = infected_only, exposures = c("non-user", "user")),
        report_deaths(cohorts, cutoff_date = cutoff_date, infected_only = infected_only, exposures = c("non-user", "irregular")),
        report_deaths(cohorts, cutoff_date = cutoff_date, infected_only = infected_only, exposures = c("non-user", "regular")),
        report_deaths(cohorts, cutoff_date = cutoff_date, infected_only = infected_only, exposures = c("irregular", "regular"))
    )
}




##' Generate a gt table for publication
##' @title overview_table
##' @param outcomes the result of overview_outcomes
##' @param cohorts a list of cohorts
##' @param cutoff_date the date after which evens should not be considered
##' @return a gt object
##' @author Greg Tucker-Kellogg

overview_gt <- function(outcome_df) {
    gc()
    gc()

    outcome_df |>
        mutate(outcome = sprintf("**%s**", stringr::str_to_title(outcome))) |>
        filter(str_ends(comparison, exposure)) |>
        select(outcome, comparison, matches("estimate|upper|lower|p_value")) |>
        mutate(across(matches("estimate"), ~ 1 - .x, .names = "{.col}_reduction")) |>
        gt(groupname_col = "outcome", rowname_col = "comparison", process_md = TRUE) |>
        tab_stubhead(label = "outcome") |>
        tab_spanner(columns = contains("simulated"), id = "A", label = md("**Simulation**")) |>
        tab_spanner(columns = contains("kerr"), id = "B", label = md("**KB22/KC22**")) |>
        tab_spanner(columns = contains("corrected"), id = "C", label = md("**'Corrected'**")) |>
        cols_move("kerr_estimate", "simulated_p_value") |>
        cols_move("kerr_upper", "kerr_estimate") |>
        cols_move("kerr_lower", "kerr_upper") |>
        cols_move("kerr_estimate_reduction", "kerr_lower") |>
        cols_move("kerr_p_value", "kerr_estimate_reduction") |>
        cols_move("simulated_estimate_reduction", "simulated_lower") |>
        cols_move("corrected_estimate_reduction", "corrected_lower") |>
        fmt_number(rows = everything(), columns = matches("estimate|lower|upper")) |>
        fmt_percent(rows = everything(), columns = ends_with("reduction"), decimals = 0) |>
        fmt_number(rows = everything(), columns = contains("p_value"), decimals = 3) |>
        sub_small_vals(rows = everything(), columns = contains("p_value"), threshold = 0.001) |>
        cols_merge_range(
            col_begin = "simulated_lower",
            col_end = "simulated_upper"
        ) |>
           cols_merge_range(
            col_begin = "kerr_lower",
            col_end = "kerr_upper"
        ) |>
        cols_merge_range(
            col_begin = "corrected_lower",
            col_end = "corrected_upper"
        ) |>
        cols_label(
            ends_with("estimate") ~ "RR",
            ends_with("reduction") ~ "Risk Red.",
            contains("lower") ~ "95% CI",
            contains("p_value") ~ "p.val"
        ) |>
        tab_stubhead(label = "Summary statistics")
}


overview_gt_to_latex <- function(gt_obj, model_name, probabilistic, infected_only,label) {
    caption <- sprintf(
        "\\caption{%s %s %s}",
        sprintf("%s model with %s stop on infection.", model_name, ifelse(probabilistic, "probabilistic", "uniform")),
        sprintf("Stop probability: %s.", ifelse(probabilistic, sprintf("%3.2f (irregular), %3.2f (regular)", max(model$p_inf_stop), min(model$p_inf_stop)), "1.0 (all individuals)")),
        sprintf(
            "Statistics for hospitalisations and deaths are %s",
            ifelse(infected_only, "limited to infected individuals.", "reported for all individuals in the cohort.")
        )
    )

    as_latex(gt_obj) |>
        as.character() |>
        str_replace_all(c(
            "[$]" = "",
            "longtable" = "tabulary",
            "0.000" = "<0.001"
        )) |>
        str_replace(
            fixed("\\begin{tabulary}{l|rrrrrrrrrrrr}"),
            "\n\n\\begin{table}[h!tb]\n\\footnotesize\n\\begin{tabulary}{\\linewidth}{lRrrRRrrRRrrR}"
        ) |>
        str_replace(
            fixed("\\begin{tabulary}{l|rrrrrrrr}"),
            "\n\n\\begin{table}[h!tb]\n\\footnotesize\n\\begin{tabulary}{\\linewidth}{lRrrRRrrR}"
        ) |>
        paste0(sprintf("%s\n%s\n\\end{table}\n\n", caption,glue("\\label{{tab:{label}}}")), collapse = "\n")
}


### Start

overview_dfs <- tryCatch( {
    read_rds(here::here("results", "overview_dfs.rds"))
},
error = function(cond) {
    message(paste("could not load file"))
    message("Here's the original error message:")
    message(cond)
    message()
    list()
}
)

model_names <- tibble(model=c('sim1','sim2','sim3'),
       name=c('i-ENR', 'i-INF','i-KC22'))


overview_table <- expand_grid(model_names,stopping=c("probabilistic","uniform"),subset=c("infected","all")) |>
    mutate(df=paste(model,subset,stopping,sep="-"))


if (length(overview_dfs) == 0) {
    needed <- overview_table
} else {
    needed <- anti_join(overview_table,tibble(df=names(overview_dfs)))
}



needed <- mutate(needed,
                 reload=ifelse(model == lag(model),FALSE,TRUE),
                 restop=(ifelse(reload | (stopping != lag(stopping)),TRUE,FALSE))) |>
    replace_na(list(reload=TRUE,restop=TRUE))

needed


### Supplementary


if (nrow(needed) > 0) { 
    for (i in 1:nrow(needed)) {
        dataset <- as.list(needed[i,])
        message(glue("*********************** Model is {dataset$df},  **********************"))
        if (dataset[["reload"]]) {
            gc()
            message(glue("------------------------ Loading {dataset$model}"))
            cohorts <- read_rds(here("results", sprintf("%s.rds", dataset$model)))
            gc()
        }
        
        if (dataset[['restop']]) {
            message(glue("----------- reassigning regularity with stopping rule = {dataset$stopping}"))
            cohorts <- reassign_regularity(cohorts, probabilistic = (dataset$stopping == "probabilistic"))
        }
        message(glue("----------- working with {dataset$subset} samples"))
        overview_dfs[[dataset$df]] <- overview_outcomes(cohorts, infected_only = (dataset$subset == "infected"))
        write_rds(overview_dfs,file=here::here('results',"overview_dfs.rds"))
    }
}

supplemental_tables <- map(setNames(model_names$name,model_names$model),
                           function(s) sprintf(
                                           "\n\n\\clearpage\n\\subsection{%s model simulations}\n\n",
                                           s))

supplemental_tables <- c(supplemental_tables,
                         map(setNames(1:nrow(overview_table),overview_table$df),
                             function(i) {
                                 l <- as.list(overview_table[i,])
                                 overview_gt(overview_dfs[[l$df]]) |>
                                     cols_hide(columns=matches("corrected")) |>
                                     overview_gt_to_latex(l$name,
                                                          probabilistic = l$stopping == "probabilistic",
                                                          infected_only = l$subset == "infected",
                                                          label=l$name)
                             }))


supplemental_tables <- supplemental_tables[sort(names(supplemental_tables))]
cat(as.character(supplemental_tables), file = here("results/supplementary.tex"))


table4 <- overview_gt(overview_dfs[["sim3-all-probabilistic"]])
cat(as.character(overview_gt_to_latex(table4,"i-KC22",
                                      probabilistic=TRUE,
                                      infected_only=FALSE,
                                      label="main_results")),
    file=here("results/table4.tex"))

   
