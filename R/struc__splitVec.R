#' @title Splitting a vector into a list of n vectors
#' @description This function splits a vector into a list of n vectors
#' @param v Input vector to be broken
#' @param n Number of vectors to be generated from v
#' @return An unnamed list of n vectors
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 3.0
# Publication: 10 June 2022; last update: 10 June 2022

splitVec <- function(v, n) {
    N <- length(v)
    if (n > 0 && n < N) {
        l <- split(x = v, cut(x = seq_along(v), breaks = n, labels = FALSE))  # Reference: stackoverflow.com/questions/3318333/split-a-vector-into-chunks
    } else {
        l <- v
    }
    return(l)
}
