#' @title Generating a list for construction a network using the function forceNetwork from the networkD3 package
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Editions: 2015, 17/3/2016

buildNetwork <- function(df, w) {
    # initialise variables
    nodes <- sort(union(df[, 1], df[, 2]))  # get a sorted list of all node names
    ID <- seq(from = 0, to = length(nodes) - 1, by = 1)
    group <- rep(1, times = length(ID))

    n <- nrow(df)
    node1 <- node2 <- val <- numeric(n)  # assign a vector of numerics
    links <- data.frame(node1, node2, val)

    # convert the dataframe into a network file
    for (i in 1 : n) {
        links$node1[i] <- as.numeric(which(nodes == df[i, 1]) - 1)  # As the ID starting from zero
        links$node2[i] <- as.numeric(which(nodes == df[i, 2]) - 1)
        links$val[i] <- df[i, 3] * w
    }

    r <- list(links, data.frame(ID, name = nodes, group))
    names(r) <- c("links", "nodes")
    return(r)
}
