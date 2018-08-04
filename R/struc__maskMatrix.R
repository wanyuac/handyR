#' @title Generating a matrix consisting of 0/1 for filtering values in a target matrix.
#' @description For example, m <- m * maskMatrix(m, items_to_keep.dataframe, keep.diag = True) will replace unwanted cells with zero.
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 8 Nov 2015

########## FILTERRING OF ELEMENTS IN A MATRIX ##########
maskMatrix <- function(m, df, default = 0, keep.diag = TRUE) {
    # df: $node1, $node2
    # m must be a synmetric matrix with row names and column names defined
    # keep.diag: whether to keep diagonal elements of the target matrix

    n <- nrow(m)
    mask <- matrix(default, nrow = n, ncol = n)
    items <- colnames(m)
    rownames(mask) <- colnames(mask) <- items
    names(df)[1 : 2] <- c("item1", "item2")

    for (i in 1 : nrow(df)) {
        a <- df$item1[i]
        b <- df$item2[i]
        mask[a, b] <- mask[b, a] <- 1
    }

    if (keep.diag) {
        for (i in 1 : nrow(mask)) {
            mask[i, i] <- 1
        }
    }

    return(mask)
}