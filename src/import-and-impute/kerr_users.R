## import and clean up the Kerr infected data

library(here)
library(tidyverse)
library(lubridate)
library(readxl)
source(here("src","data","kerr_datasets.R"),local=TRUE)

kerr_users <- read_excel(kerr_files$ds1) %>%
    janitor::clean_names() %>%
    mutate(date_birth=as.Date(birth_date),exposure='user') %>%
    select(-birth_date) %>% 
    arrange(date_birth) 

