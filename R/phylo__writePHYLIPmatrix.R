#' @title Create a PHYLIP-formatted distance matrix from the output of function makeFastANIDistMatrix
#' @description This function takes as input a square matrix and write a PHYLIP-formatted matrix file.
#' See http://evolution.genetics.washington.edu/phylip/doc/distance.html for an example of the format.
#' @param M The distance matrix with row and column names
#' @param f Path of the output file
#' @return NULL
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# Copyright 2020 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 3.0
# Publication: 6 April 2020

writePHYLIPmatrix <- function(M, f) {
    n <- nrow(M)  # Number of genomes
    ids <- rownames(M)
    write(paste0("\t", n), file = f)  # The write function appends a newline character automatically.
    for (i in 1 : n) {
        write(paste(c(ids[i], as.character(M[i, ])), collapse = "\t"), file = f, append = TRUE)
    }
}
