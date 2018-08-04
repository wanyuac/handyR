#' @title Making a contingency table from four entries
#' @description This programme stimulates a two-by-two contingency table.
#' Structure of the contingency table
#'       B
#'       1  0
#' A  1  a  b
#'    0  c  d
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 21 Oct 2015

#========== CONSTRUCT A CONTINGENCY TABLE ==========
mkContingencyTable <- function(a, b, c, d) {
  m <- matrix(data = c(a, c, b, d), nrow = 2, ncol = 2)
  rownames(m) <- colnames(m) <- c("presence", "absence")
  return(m)
}

#========= FISHER'S EXACT TEST ==========
#fisher.test(mkContingencyTable(a = 1, b = 100, c = 100, d = 399), alternative = "two.sided")