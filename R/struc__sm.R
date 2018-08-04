#' @title Show a diagonal matrix m[start : end, start : end].
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 14 Dec 2016

sm <- function(m, start = NULL, end) {  # show matrix
	if (is.null(start)){
		start <- 1
	}
    print(m[start : end, start : end])
}
