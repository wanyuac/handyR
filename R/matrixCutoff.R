#' @title Resets a value in an input matrix if it is smaller/greater than the threshold.
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licence: Apache License, Version 2.0
# Development history: 27 Sept 2015

matrixCutoff <- function(m, cutoff, comp = ">", new.val = 0) {
  for (i in 1 : nrow(m)) {
    for (j in 1 : ncol(m)) {
      if (comp == ">") {
        if(m[i, j] > cutoff) {
          m[i, j] <- new.val
        }
      } else {
        if(m[i, j] < cutoff) {
          m[i, j] <- new.val
        }
      }
    }
  }
  return(m)
}