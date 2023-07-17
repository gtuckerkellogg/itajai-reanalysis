## Return a tibble of Kerr et al Data Set 2, the adult participants with usage, but cleaning up invalid usage
## entries

library(here)

kerr_adult_participants <- local({
    source(here("src","data","global_params.R"),local=TRUE)
    source(here("src","data","kerr_datasets.R"),local=TRUE)
                                        # maximum # of pills is 4/day. Some of the raw entries are absurd
    max_pills <- 4*length(c(seq(1,paper$duration,by=15),seq(2,paper$duration,by=15)))

    d <- readxl::read_excel(kerr_files$ds2) |>
        janitor::clean_names() |>
        dplyr::rename(total_tablets='accumulated_number_of_tablets_each_tablet_6mg')
    invalid_rows <- which(d$total_tablets > 80)
    d$total_tablets[invalid_rows] <- sample(d$total_tablets[-invalid_rows],
                                            length(invalid_rows),replace=TRUE)
    d
})




