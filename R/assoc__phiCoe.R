#' @title Calculating a phi association coefficient from two binary variables.
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 8 Nov 2015

phiCoe <- function(x, y = NULL)
{
  # This function was copied from the package GenomicRanges, which is not supporting R/3.1.1.
  if (is.null(y)) {
    if (!is.integer(x) || length(x) != 4L)
      stop("when 'y' is not supplied, 'x' must be ",
           "a 2x2 integer matrix filled by columns or an integer vector of length 4")
    # The contingency table is filled by columns.
    a <- x[1L]
    c <- x[2L]
    b <- x[3L]
    d <- x[4L]
  } else {
    if (!is.logical(x) || !is.logical(y) || length(x) != length(y))
      stop("when 'y' is supplied, 'x' and 'y' must be ",
           "2 logical vectors of the same length")
    a <- sum(x  &  y)
    b <- sum(x  & !y)
    c <- sum(!x &  y)
    d <- sum(!x & !y)
  }
  a <- as.double(a)  # to avoid the problem of "integer overflow"
  b <- as.double(b)
  c <- as.double(c)
  d <- as.double(d)
  div <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
  return((a * d - b * c) / div)
}
