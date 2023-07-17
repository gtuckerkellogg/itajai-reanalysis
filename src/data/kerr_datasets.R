library(here)

kerr_files <- purrr::map(1:7,function(n) { Sys.glob(file.path(here("data/kerr"), sprintf("Dataset %d*.xlsx",n))) })

names(kerr_files) <- purrr::map_chr(1:7, function(n) sprintf("ds%d",n))
