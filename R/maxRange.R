#' @title Obtain the maximum range comprising two numeric vectors v1 and v2
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 16 Nov 2015

maxRange <- function(v1, v2) {
	r1 <- range(v1)
	r2 <- range(v2)
  lower.bound <- min(r1[1], r2[1])
  upper.bound <- max(r1[2], r2[2])
  return(c(lower.bound, upper.bound))
}