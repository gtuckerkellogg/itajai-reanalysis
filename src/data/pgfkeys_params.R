#!/usr/bin/env Rscript

local({

    library(lubridate)
    library(purrr)
    library(here)
    library(stringr)
    source(here::here("src/data/global_params.R"), local = TRUE)

    keysfile <- here::here("results","keys.tex")
    invisible(suppressWarnings(file.remove(keysfile)))
    keysfile <- file(keysfile,"w")

    cat(
        map_chr(sort(names(as.list(paper))),
                function(name) {
                    sprintf("\\pgfkeyssetvalue{/KC22/%s}{%s}\n",name,paper[[name]])
        }
        ),
        file=keysfile)

    cat(c(sprintf("\\pgfkeyssetvalue{/model/%s}{%s}\n","high_stop",min(model$p_inf_stop)),
          sprintf("\\pgfkeyssetvalue{/model/%s}{%s}\n","low_stop",max(model$p_inf_stop))),
        file=keysfile)

    if (!exists('results')) {
        result=list() 
    }
    
    cat(
        map_chr(sort(names(as.list(results))),
                function(name) {
                    val = results[[name]]
                    sprintf("\\pgfkeyssetvalue{/results/%s}{%s}\n",name,val)
        }
        ),
        file=keysfile)

    cat(paste("\n",
              "\\newcommand{\\themodel}[1]{\\pgfkeysvalueof{/model/#1}}",
              "\\newcommand{\\result}[1]{\\pgfkeysvalueof{/results/#1}}",
              "\\newcommand{\\numresult}[1]{\\num{\\pgfkeysvalueof{/results/#1}}}",
              "\\newcommand{\\kerr}[1]{\\pgfkeysvalueof{/KC22/#1}}",
              "\\newcommand{\\numkerr}[1]{\\num{\\pgfkeysvalueof{/KC22/#1}}}",
              "\n",
              sep="\n"),
        file=keysfile)

    rm(keysfile)
})
