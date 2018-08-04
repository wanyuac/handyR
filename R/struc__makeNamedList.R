#' @title Returns a one- or two-dimensional list with names for each dimension.
#' @param names.1 The first vector of characters as names for each level of the list
#' @param names.2: The second vector of characters as names for each level of the list.
#' e.g., an element in such a list can be assessed using x[["L1"]][["L2"]].
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 16/3/2016, 28/4/2016

makeNamedList <- function(names.1 = NA, names.2 = NA) {
    # names.1 and names.2: names for the first and second levels of elements in a list
    # Each element in a two-level list can be accessed by x[[names.1]][[names.2]]
    if (!is.na(names.1)[1]) {
        if (!is.na(names.2)[1]) {
            # create a two-level list if both levels of names are present
            x <- vector(mode = "list", length = length(names.2))
            names(x) <- names.2
            y <- vector(mode = "list", length = length(names.1))
            names(y) <- names.1
            for (item in names.1) {
                y[[item]] <- x  # y[[names.1]][[names.2]]
            }
        } else {
            # create a usual list if names for the second-level elements are not configured
            y <- vector(mode = "list", length = length(names.1))
            names(y) <- names.1
        }
    } else {
        if (!is.na(names.2)[1]) {
            # the same as the situation where only names.1 is configured
            y <- vector(mode = "list", length = length(names.2))
            names(y) <- names.2
        } else {
            # A NULL value is returned if neither level is configured.
            y <- NULL  # no output
        }
    }
    return(y)
}

