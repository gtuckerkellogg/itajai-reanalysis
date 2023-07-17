## Provide a cleaned up dat_all and citywide_infected. The citywide_infected is filtered to Itajai,
## detailed in a comment below

library(here) #independent of Rstudio, sets here() to the git root
library(tidyverse)
library(readr)


dat_all <- local({
    ## Now cross-check with Brazilian state data
    ## download this data from google drive, too big for git
    INFLUD20 <- here('data',"INFLUD20-26-09-2022.csv")
    latest_birth_date <- as.Date("2020-01-01") - years(18)

    ## The data is a bit messy, should clean it up 
    colspec <- spec_delim(INFLUD20,delim=";")

    for (i in which(str_starts(names(colspec$cols),"DT_"))) {
        colspec$cols[[i]] <- col_date(format="%d/%m/%Y")
    }

    colspec$cols[['OUT_ANIM']] <- col_character()
    colspec$cols[['FLUASU_OUT']] <- col_character()
    colspec$cols[['FLUBLI_OUT']] <- col_character()

    dat_all <- suppressWarnings(read_delim(INFLUD20,
                                           delim=';',col_types=colspec))
                                        # print(problems(dat_all),n=30) ## seems OK now

    dat_all <- dat_all[-problems(dat_all)$row,] |>
        mutate(hospitalised=HOSPITAL==1,
               death=EVOLUCAO>1) |>
        rename(date_birth=DT_NASC,
               date_onset=DT_SIN_PRI,
               date_notified=DT_NOTIFIC,
               date_hospitalised=DT_INTERNA,
               date_death=DT_EVOLUCA,
               sex=CS_SEXO) |>
        replace_na(list(hospitalised=FALSE,death=FALSE)) |>
        filter(date_birth <= latest_birth_date,
               !(death & is.na(date_death)),
               !(hospitalised & is.na(date_hospitalised)))
    dat_all
})


citywide_infected <- local({

## How to filter.  Email from Ana, 1 Dec 2022:
## ID_MUNICIP : the city where the case has been notified to the National Health Authority.

## ID_REGIONA: micro-region where the case has been notified to the National Health Authority --> basically
## means Itajaí city and other surrounding small towns.

## ID_MN_INTE  and ID_RG_INTE : respectively, the city and micro-region where the medical center that admitted
## the patient was located.


## ID_MN_RESI and ID_RG_RESI : respectively, the city and micro-region where the patient LIVED


## Kerr and cols. could neither retrieve data from medical centers outside Itajaí nor included patients from
## other cities even if they had taken ivermectin from the Itajaí's citywide program. Therefore, I would keep
## the last 4 filters and ignore "ID_MUNICIP" and/or "ID_REGIONA" since the case notification could come from
## the first hospital visit in a center outside of Itajaí city - the patients would be eventually transferred
## to Itajaí city.

                                        # Ana's recommendation
dat_all |>
    filter(PCR_SARS2 == 1 | CLASSI_FIN==5) |>    
    filter_at(vars(any_of(c("ID_MN_INTE","ID_RG_INTE","ID_MN_RESI","ID_RG_RESI"))),
              any_vars( . == "ITAJAI"))

})
