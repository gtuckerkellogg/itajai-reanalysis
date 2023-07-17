library(here)

kerr_files <- list(ds1="Dataset 1 - All participants of the Itajai program (n=147,223).xlsx",
             ds2="Dataset 2 - Participants above 18 yo not infected before the program (Before July 7, 2020) from the city of Itajaí (n=113,844).xlsx",
             ds3="Dataset 3 - Participants and non-p articipants above 18 yo infected during the period of the program from Itajai and outside (n=7,345).xlsx",
             ds4="Dataset 4 - Participants above 18 yo infected during  the program from Itajaí and outside (n= 4,311).xlsx",
             ds5="Dataset 5 - Participants above 18 yo infected during the program from Itajaí only (n=4,194).xlsx",
             ds6="Dataset 6 - Participants above 18 years old infected during the program from outside Itajaí (n=117).xlsx",
             ds7="Dataset 7 - Non-participants above 18 yo infected during the period of the program (n=3,034).xlsx") |>
    purrr::map(~here("data/kerr",.x))





