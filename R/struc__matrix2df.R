#' @title Coverting the upper triangle of a symmetric numeric matrix into a data frame of three variables.
#' @description For speed, this function does not check symetry of the input matrix. It is a user's responsibility
#' to guarantee the input is a symmetric matrix. In addition, row names and column names must be
#' the same to get a correct output.
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 27 Sept 2015, 29/4/2016, 12/7/2017 (substantially optimised the algorithm to make it dramatically faster than before under the "small" mode)

matrix2df <- function(m, diag = FALSE, named = FALSE, small = TRUE, val.class = "numeric") {
    # m: a symmetric matrix of integers;
    # diag: whether the diagnoal should be included.
    # named: whether the matrix is named on both columns and rows.
    # small: whether the input matrix is small in dimenstions
    # val.class: class (integer, numeric or character) of values

    if (named) {
        ids <- colnames(m)  # assuming that row names = column names
        init.id.class <- character
        if (sum(rownames(m) != ids) > 0) {
            print("Warning: rearranging rows to match column names.")
            m <- m[ids, ]
        }
    } else {
        ids <- 1 : ncol(m)  # use indices instead
        init.id.class <- integer
    }

    if (val.class == "numeric") {
        conv.class <- as.numeric
        init.class <- numeric
    } else if (val.class == "integer") {
        conv.class <- as.integer
        init.class <- integer
    } else {
        conv.class <- as.character
        init.class <- character
    }

    n <- ncol(m)  # Note that n == NULL if m is a single numeric, then character(NULL) returns an error of invalid "length" argument.
    if (diag) {
        row.start <- 1
        col.start <- 1
        row.end <- n
        col.end <- n
    } else {
        row.start <- 1
        col.start <- 2
        row.end <- n - 1
        col.end <- n
    }

    if (small) {  # The small-matrix mode: row names and column names will go into the data frame
        df <- data.frame(V1 = init.id.class(0), V2 = init.id.class(0), v = init.class(0), stringsAsFactors = FALSE)

        # turn off the "small" option when m is extremely large
        ids.row <- ids[row.start : row.end]  # n or n - 1 elements
        for (i in ids.row) {
            col.indices <- col.start : col.end  # column indices to be read from a row
            ids.col <- ids[col.indices]  # numerics or characters
            r <- conv.class(m[i, col.indices])  # read length(ids.col) values from a row
            df <- rbind.data.frame(df, data.frame(V1 = rep(i, times = length(ids.col)),
                                                  V2 = ids.col, v = r, stringsAsFactors = FALSE),
                                   stringsAsFactors = FALSE)
            col.start <- col.start + 1
        }
    } else {
        # The large-matrix mode, which is a memory-efficient way for processing an extremely large matrix.
        # A list of two variables is returned in this mode instead of a data frame.
        # Values in the matrix are transferred to a single vector, where row names or column names are discarded.
        # The second vector is included in the list in order to help users to access row names and column names.
        # The second vector records the index of a value in the first vector, where a new row in the matrix starts being read.
        k <- 1  # initialise a counter for values transferred
        N <- (row.end + row.start) * (row.end - row.start + 1) / 2  # number of values to be transferred
        v <- integer(N)  # make a vector of integers for values to be read
        row.head <- integer(row.end - row.start + 1)  # An index is recorded for the start of each row. So the vector's length equals the number of rows to be read.
        for (i in row.start : row.end) {
            row.head[i] <- k
            for (j in col.start : col.end) {
                v[k] <- m[i, j]
                k <- k + 1
            }
            col.start <- col.start + 1
        }
        df <- list(v = v, row.head = row.head)
    }

    return(df)
}
