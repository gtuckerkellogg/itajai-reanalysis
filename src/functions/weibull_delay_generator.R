## Given an objective mean and sd, return a discrete delay generator function based on a
## corresponding Weibull distribution. The Weibull shape and scale parameters are
## optimised using the provided mean and sd.

library(rlang)
library(here)

delay_generator <- local({
    source(here("src/data/global_params.R"),local=TRUE)        
    function(.data,dt1,dt2,min_diff=1,max_diff=model$delay_limit) {

        time_difference <- na.omit(as.numeric(.data[[dt2]] - .data[[dt1]]))
        time_difference <- time_difference[between(time_difference,min_diff,max_diff)]
        objective_mean <- mean(time_difference,na.rm=TRUE)
        objective_sd <- sd(time_difference,na.rm=TRUE)
        objective <- function (weibull_params) {
            scale = weibull_params[1]
            shape = weibull_params[2]
            weibull_mean = scale * gamma(1 + 1/shape)
            weibull_sd = sqrt((scale^2) * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2))
            return((objective_mean - weibull_mean)^2 + (objective_sd - weibull_sd)^2)
        }
        opt <- optim(c(1,1),objective)
        scale <- opt$par[1]
        shape <- opt$par[2]
        list(scale=opt$par[1],shape=opt$par[2])
        function(n) { round(rweibull(n,shape=shape,scale=scale)) }
    }
    })
