#' @title Make a distance matrix from FastANI output (tab-delimited file)
#' @description This function coverts a FastANI output into a distance matrix.
#' @param f Path to the tab-delimited output file of FastANI
#' @param keep_asym A logical flag specifying whether to keep the original asymmetric distance matrix.
#' @param frac A logical flag specifying whether to convert percentages to decimal fractions. This option
#' does not affect the tree topology as the change is proportional.
#' @return One or two n-by-n distance matrices (depending on keep_asm), where n denotes the number of genomes.
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# Copyright 2020 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 3.0
# Publication: 6 April 2020

makeFastANIDistMatrix <- function(f, keep_asym = FALSE, frac = FALSE) {
    # Initiation
    ani <- read.delim(file = f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1 : 3]
    names(ani) <- c("Query", "Reference", "ANI")
    ids <- sort(union(ani$Query, ani$Reference), decreasing = FALSE)
    ani$D <- 100 - ani$ANI  # Calculate distances from ANIs
    ani <- ani[, -3]  # Remove the column "ANI"

    if (frac) {
        ani$D <- ani$D / 100  # Convert percentages to decimal fractions
        precision <- 6  # Number of decimals to keep
    } else {
        precision <- 4  # The same as FastANI
    }

    n <- length(ids)
    M <- matrix(data = NA, nrow = n, ncol = n, dimnames = list(ids, ids))
    diag(M) <- 0

    # Stage one: copy values from the data frame to matrix M
    for (i in 1 : nrow(ani)) {
        rw <- ani[i, ]  # Extract one row from the data frame to increase the speed
        q <- rw$Query
        r <- rw$Reference
        if (r != q) {
            M[q, r] <- rw$D
        }
    }

    # Stage two: convert M into a symmetric matrix by taking the mean distance between every pair of genomes
    # This is the same method that FastANI uses for generating the PHYLIP-formatted lower triangular matrix.
    # See https://github.com/ParBLiSS/FastANI/issues/36
    if (keep_asym) {
        M_asym <- M
    }

    for (i in 1 : (n - 1)) {
        for (j in (i + 1) : n) {
            val_up <- M[i, j]  # The value in the upper triangle
            val_lo <- M[j, i]  # The value in the lower triangle
            v <- round((val_up + val_lo) / 2, digits = precision)  # The same precision as FastANI (after dividing values by 100)
            M[i, j] <- v
            M[j, i] <- v
        }
    }

    # Return the result
    if (keep_asym) {
        out <- list("D" = M, "D_asym" = M_asym)
    } else {
        out <- M
    }

    return(out)
}
