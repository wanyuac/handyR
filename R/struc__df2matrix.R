#' @title Coverts a dataframe that corresponds to either the upper/lower triangle in a matrix into a matrix.
#' @description Variables in the data frame: c(name1, name2, value)
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 10 Oct 2015, 12/2/2017

df2matrix <- function(df, diag = 1, dimension = NULL, symmetric = FALSE, replace.na = NA) {
    # diag: default values of diagonal cells;
    # replace.na: alternative values for matrix cells that are not covered by the data frame
    # Variables in the input data frame: node1  node2  value
    # symmetric: is the target matrix symmetric?
    # dimension: c(nrow, ncol) for a user-specified dimension of the target matrix

    # initialises a matrix
    if (length(dimension) == 2) {
        m <- matrix(data = replace.na, nrow = dimension[1], ncol = dimension[2])
        rownames(m) <- as.character(1 : nrow(m))
        colnames(m) <- as.character(1 : ncol(m))
        n <- min(dimension)
    } else {
        elements <- sort(union(df[ , 1], df[ , 2]))  # colnum names and row names
        n <- length(elements)
        dimension <- c(n, n)
        m <- matrix(data = replace.na, nrow = n, ncol = n)
        rownames(m) <- colnames(m) <- elements
    }

    # fills the diagonal
    for (i in 1 : n) {
        m[i, i] <- diag
    }
    filled <- n

    # fills cells with values recorded in the input data frame
    if (symmetric) {
        for (i in 1 : nrow(df)) {
            node1 <- df[i, 1]
            node2 <- df[i, 2]
            m[node1, node2] <- m[node2, node1] <- df[i, 3]
            filled <- filled + 2
        }
    } else {
        for (i in 1 : nrow(df)) {
            node1 <- df[i, 1]
            node2 <- df[i, 2]
            m[node1, node2] <- df[i, 3]
            filled <- filled + 1
        }
    }

    # check if any values remain missing
    #if (sum(is.na(m)) > 0) {
    #    print("Warning: some cells have not been filled.")
    #}
    if (dimension[1] * dimension[2] > filled) {
        print("Warning: some cells have not been filled.")
    }

    # fills other cells not described in the input data frame
    #if (!is.na(replace.na)) {
    #    cat("Replacing NAs with ", replace.na, ".\n", sep = "")
    #   for (i in 1 : n) {
    #        for (j in 1 : n) {
    #            if(is.na(m[i, j])) {
    #                m[i, j] <- m[j, i] <- replace.na
    #            }
    #        }
    #    }
    #}

    return(m)
}