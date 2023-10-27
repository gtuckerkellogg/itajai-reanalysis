## Provide a cleaned up dat_all and citywide_infected. The citywide_infected is filtered to Itajai,
## detailed in a comment below

library(here) #independent of Rstudio, sets here() to the git root
library(tidyverse)
library(readr)
source(here("src/data/global_params.R"), local = TRUE)
source(here("src/data/keep_results.R"), local = TRUE)



if (!exists("dat_all")) { 

    dat_all <- local({

        race_dict <- c("1" = "Branca",
                       "2" = "Preta",
                       "3" = "Amarela",
                       "4" = "Parda",
                       "9" = "Ignorado")
        
        dat_all_file <- here("data","dat_all.rds")
        
        if (!file.exists(dat_all_file)) { 
            ## Now cross-check with Brazilian state data
            ## download this data from google drive, too big for git
            INFLUD20 <- here('data',"INFLUD20-26-09-2022.csv")
            ## Optional newer file, misses a couple of matches
                                        #  INFLUD20 <- here('data',"INFLUD20-01-05-2023.csv")
            latest_birth_date <- as.Date("2020-01-01") - years(18)

            ## The data is a bit messy, should clean it up 
            colspec <- spec_delim(INFLUD20,delim=";")

            for (i in which(str_starts(names(colspec$cols),"DT_"))) {
                colspec$cols[[i]] <- col_date(format="%d/%m/%Y")
            }

            colspec$cols[['OUT_ANIM']] <- col_character()
            colspec$cols[['CS_RACA']] <- col_character()
            colspec$cols[['FLUASU_OUT']] <- col_character()
            colspec$cols[['FLUBLI_OUT']] <- col_character()

            dat_all <- suppressWarnings(read_delim(INFLUD20,
                                                   delim=';',col_types=colspec))
                                        # print(problems(dat_all),n=30) ## seems OK now

            dat_all <- dat_all[-problems(dat_all)$row,] |>
                rename(date_birth=DT_NASC,
                       date_onset=DT_SIN_PRI,
                       date_notified=DT_NOTIFIC,
                       race=CS_RACA,
                       date_hospitalised=DT_INTERNA,
                       date_death=DT_EVOLUCA,
                       type_2_diabetes=DIABETES,
                       sex=CS_SEXO,
                       hospitalised=HOSPITAL) |>
                mutate(hospitalised = ((hospitalised == 1) | (!is.na(date_hospitalised))),
                       type_2_diabetes = ((type_2_diabetes == 1)),
                       race=unname(race_dict[race]),
                       death = ((EVOLUCAO == 2))) |> 
                replace_na(list(hospitalised=FALSE,death=FALSE,type_2_diabetes=FALSE)) |>
                filter(date_birth <= latest_birth_date)

                write_rds(dat_all,dat_all_file)
        } else { 
            dat_all <- read_rds(dat_all_file)
        }
        dat_all
    })
}

if (!exists("citywide_infected")) {
    citywide_infected <- local({
        
        ## We consider first the possible matches from DATASUS, being as generous as possible considering the
        ## wording of the paper.
        
        ## The city-wide program was designed for residents of \Itajai, which we capture by filtering on
        ## resident status.
        
        ## For non-residents who joined the program, we also include ID_MUNICIP and IT_MN_INTE, assuming that
        ## they would have been recorded through the local hospitals or reported by the city itself.
        
        citywide_infected <-
            dat_all |>
            filter(PCR_SARS2 == 1 | CLASSI_FIN==5) |>
            filter(date_onset >= as.Date("2020-01-01"),date_onset <= as.Date("2020-12-31")) |>
            filter(if_all(c("ID_RG_RESI","ID_MN_RESI"), ~ . == "ITAJAI") |    # True city residents
                   if_any(c("ID_MUNICIP","ID_MN_INTE"),~ . == "ITAJAI")) # Reported through city

                   citywide_infected

    })
}
