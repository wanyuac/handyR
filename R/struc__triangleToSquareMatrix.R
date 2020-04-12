#' @title Converting a data frame representing a triangle matrix into a square matrix.
#' @description This function retrieves IDs from the first two columns and values from the third column.
#' Note that this function assumes the diagonal of the output matrix is comprised of a constant.
#' @param d Input data frame
#' @param diag_val The value on the diagonal (Default: NA)
#' @param sort_id A logical flag specifying whether to sort row and column names of the output matrix
#' @param sort_decreasing A logical flag specifying whether to sort row and column names in a decreasing order
#' @return A synmetric matrix
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# Copyright 2020 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 3.0
# Publication: 12 April 2020

triangleToSquareMatrix <- function(d, diag_val = NA, sort_id = FALSE, sort_decreasing = FALSE) {
    # Initialisation
    if (ncol(d) > 3) {
        d <- d[, 1 : 3]
    }

    names(d) <- c("x", "y", "v")
    ids <- union(d$x, d$y)
    n <- length(ids)

    if (sort_id) {
        ids <- sort(ids, decreasing = sort_decreasing)
    }

    M <- matrix(NA, nrow = n, ncol = n, dimnames = list(ids, ids))

    if (!is.na(diag_val)) {
        diag(M) <- diag_val
    }

    # Populate the matrix
    for (i in 1 : nrow(d)) {
        r <- d[i, ]
        x <- d$x
        y <- d$y
        v <- d$v
        M[x, y] <- v
        M[y, x] <- v
    }

    return(M)
}