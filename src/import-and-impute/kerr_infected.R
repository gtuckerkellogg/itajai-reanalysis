## import and clean up the Kerr infected data

library(here)
library(tidyverse)
library(lubridate)


kerr_infected <- local({
    source(here("src","data","kerr_datasets.R"),local=TRUE)
    readxl::read_excel(kerr_files$ds3) |>
        janitor::clean_names() %>%
        mutate(date_birth=as.Date(birth_date)) %>%
        select(-birth_date) %>% 
        rename(hospitalised=matches("hospitalization"),
               death=matches('death'),
               race=matches('race')) %>%
        mutate(hospitalised=hospitalised==2,
               death=death==1,
               type_2_diabetes=as.logical(type_2_diabetes)) %>%
        arrange(date_birth)
})



