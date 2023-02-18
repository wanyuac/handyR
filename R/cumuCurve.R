#' @title Calculating cumulative proportions of objects in a sample population
#' @description This function calculates cumulative proportions of objects in a sample population.
#' There are two input numeric vectors, and both vectors are binned using the same width (for example, binned by 1000).
#' Therefore, all inputs are counts.
#' @param obj The object whose cumulative curve is of our interest.
#' @param pop The population from which individual objects are sampled.
#' Initially, this function is designed for processing counts in the results of the function hist(...). In this scenario,
#' both obj and pop are counts.
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018-2023 Yu Wan
# Licensed under the Apache License, Version 2.0
# First version: 19/4/2016, latest update: 18/2/2023

cumuCurve <- function(obj, pop) {
    n <- length(obj)
    cumu <- data.frame(obj = integer(n), total = integer(n), prop = numeric(n))
    if (length(pop) != n) {
        stop("Failure in pairing data: interested objects and the population have different bin numbers.")
        geterrmessage()
    } else {
        n.obj <- n.pop <- 0  # initialise variables
        # calculate cumulative proportions
        for (i in 1 : n) {
            n.obj <- n.obj + obj[i]
            n.pop <- n.pop + pop[i]
            cumu$obj[i] <- n.obj
            cumu$total[i] <- n.pop
            cumu$prop[i] <- n.obj / n.pop
        }
    }
    return(cumu)
}
