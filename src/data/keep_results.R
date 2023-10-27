library(here)

.resultsfile <- here("results/results.RData")

if (!exists('results')) {
    results <- tryCatch({
        load(.resultsfile)
        results
    },
    error = function(cond) {
        message(paste("could not load file"))
        message("Here's the original error message:")
        message(cond)
        message()
        new.env()
    }
    )

    save(results, file = .resultsfile)
    reg.finalizer(
        results,
        function(e) {
            save(results, file = .resultsfile)
        },
        onexit = TRUE)
}

.Last <- function() {
    source(here("src/data/pgfkeys_params.R"),local=TRUE)

}
