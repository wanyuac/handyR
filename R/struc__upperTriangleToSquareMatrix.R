#' @title Converting a data frame representing an upper triangle matrix (without diagonal entries)
#' into the original square matrix.
#' @description This function retrieves IDs from the first two columns and values from the third
#' column. Note that this function assumes the diagonal of the output matrix is comprised of a
#' constant.creases.
#' @param d Input data frame
#' @param diag_val The value on the diagonal (Default: NA)
#' @return A synmetric matrix
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# Copyright 2020 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 3.0
# Publication: 12 April 2020; last modification: 13 April 2020

upperTriangleToSquareMatrix <- function(d, diag_val = 0) {
    # The algorithm implemented in this function is much faster than the naive M[x, y] <- v; M[y, x] <- v method.
    # This algorithm fills a row and a column almost simultaneously.

    # Initialisation #########################
    if (ncol(d) > 3) {
        d <- d[, 1 : 3]
    }
    names(d) <- c("x", "y", "v")

    # Get row and column IDs sorted according to the occurrence count
    # This step is crucial for the efficiency of this function.
    ids <- union(d$x, d$y)
    n <- length(ids)
    if (nrow(d) != ((n ^ 2 - n) / 2)) {
        stop("Error: row number of the input data frame does not form an upper triangle matrix.")
    }

    id_nth <- setdiff(x = ids, y = unique(d$x))  # The missing ID is the n-th ID, which defines M[n, n] and is not recorded in the D.
    id_freq <- table(d$x)  # Get occurence counts of IDs
    id_freq <- data.frame(ID = names(id_freq), Freq = as.integer(id_freq), stringsAsFactors = FALSE)
    id_freq <- id_freq[order(id_freq$Freq, decreasing = TRUE), ]
    ids <- append(id_freq$ID, id_nth)  # Recover the ID order of the original upper triangle matrix.
    if (n > length(ids)) {
        stop("Error: occurrence counts of IDs are incorrect.")
    }
    M <- matrix(NA, nrow = n, ncol = n, dimnames = list(ids, ids))
    diag(M) <- diag_val

    # Populate the matrix #########################
    # Do not use M[x, y] = v; M[y, x] = v because this algorithm is extremely time-consuming
    # when n is large. Instead, the function populates the matrix by rows.
    fill_start <- 2
    for (r in id_freq$ID) {
        mark <- d$x == r  # Mark all rows of d that x == r
        ds <- d[mark, ]  # The data frame consisting of selected rows
        d <- d[!mark, ]  # The remaining data frame. It is shrinking in order to save time.
        ds$Freq <- id_freq$Freq[match(x = ds$y, table = id_freq$ID)]
        ds <- ds[order(ds$Freq, decreasing = TRUE), ]  # Sort ds in the order of the original square matrix
        v_cols <- ds$v  # Values of M[r, ] on the right of the diagonal
        fill_indices <- fill_start : n
        M[r, fill_indices] <- v_cols
        M[fill_indices, r] <- v_cols
        fill_start <- fill_start + 1
    }

    return(M)
}
