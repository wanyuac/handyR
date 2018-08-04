#' @title Taking a subset of a data frame in accordance with one or more conditions provided in the second data frame.
#' @param subject: the data frame to be subset
#' @param conditions: a data frame whose columns provide conditions for selecting rows in the subject.
#' @param col.subject: column names or indices in the subject data frame to be searched for
#' @param col.con: column names or indices in the condition data frame to be used as conditions for searching rows in the subject data frame
#'
#' Note that col.subject and col.con must be matched in terms of the order of comparisons.
#' x$a1 will be searched for values in y$a1, and so forth:
#' d <- subsetDfCon(subject = x, conditions = y, col.subject = c("a1", "a2), col.con = c("a1", "b1"))
#' x[, 2] will be searched for values in y[, 3], for example:
#' d <- subsetDfCon(subject = x, conditions = y, col.subject = c(1, 2, 3), col.con = c(1, 3, 5))
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 19/4/2016

subsetDfCon <- function(subject, conditions, col.subject = 1, col.con = 1) {
    n.sub <- length(col.subject)  # number of columns to be searched for in the subject data frame
    n.con <- length(col.con)  # number of conditions used for subsetting the subject
    d <- NA  # for the result
    if (n.sub != n.con | n.sub > ncol(subject) | n.con > ncol(conditions)) {
        stop("Invalid configuration of conditions or target columns in the subject.")
        geterrmessage()
    } else {
        selected <- rep(FALSE, times = nrow(subject))  # Initially, no rows in the subject are selected.
        for (i in 1 : nrow(conditions)) {  # go through every row of the condition data frame
            tmp <- rep(TRUE, times = nrow(subject))  # Intially, all rows in the subject are selected
            for (j in 1 : n.con) {  # go through every condition
                tmp <- tmp & subject[, col.subject[j]] == conditions[i, col.con[j]]
            }
            selected <- selected | tmp  # add these rows into selection
        }
        d <- subject[selected, ]
    }
    return(d)
}
