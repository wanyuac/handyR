#' @title Simpson's similarity coefficient between two binary vectors
#'
#' @param x The first binary vector
#' @param y The second binary vector
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2019 Yu Wan
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 15 Feb 2019

simCoe <- function(x, y) {
    x <- as.logical(x)
    y <- as.logical(y)
    return(round(sum(x & y) / min(sum(x), sum(y)), digits = 8))
}