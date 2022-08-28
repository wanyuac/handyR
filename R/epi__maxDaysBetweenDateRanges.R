#' @title Calculating maximum days between collection dates of isolates
#' @description This function works on the output of function estimateCollectionDates.
#' @param dates A data frame containing three columns: Isolate, Collection_date_precision, Collection_date_min, Collection_date_max.
#' @return A matrix of days between collection dates of each pair of isolates
# @return A list of two elements: 1. a matrix 'Tm' of days; 2. a data frame 'Days' from the upper triangle of Tm, combining with date ranges of each isolate.
#' @export maxDaysBetweenDateRanges
#' @export `%-%`
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 21 August 2022; last update: 28 August 2022

maxDaysBetweenDateRanges <- function(dates) {
    dates <- dates[, c("Isolate", "Collection_date_precision", "Collection_date_min", "Collection_date_max")]
    dates <- dates[order(dates$Isolate, decreasing = FALSE), ]
    n <- nrow(dates)
    Tm <- matrix(data = NA, nrow = n, ncol = n, dimnames = list(dates$Isolate, dates$Isolate))  # I commoned out these for speed. The function runs very slowly when the sample size becomes large.
    #Days <- data.frame(Isolate_1 = character(), Isolate_2 = character(), Pair = character(), Days_max = integer(),
    #                   Prec_1 = character(), Prec_2 = character(), Te_1 = character(), Tl_1 = character(),
    #                   Te_2 = character(), Tl_2 = character(), stringsAsFactors = FALSE)
    for (i in 1 : (n - 1)) {
        for (j in (i + 1) : n) {
            r_i <- dates[i, ]
            r_j <- dates[j, ]
            s_i <- r_i$Isolate
            s_j <- r_j$Isolate
            t <- .maxDays(te_i = r_i$Collection_date_min, tl_i = r_i$Collection_date_max,
                          te_j = r_j$Collection_date_min, tl_j = r_j$Collection_date_max,
                          p_i = r_i$Collection_date_precision, p_j = r_j$Collection_date_precision)
            Tm[s_i, s_j] <- t
            Tm[s_j, s_i] <- t
            #Days <- rbind.data.frame(Days, data.frame(Isolate_1 = s_i, Isolate_2 = s_j,
            #                                          Pair = paste(sort(c(s_i, s_j), decreasing = FALSE), collapse = "/"),
            #                                          Days_max = .maxDays(te_i, tl_i, te_j, tl_j, p_i, p_j),
            #                                          Prec_1 = p_i, Prec_2 = p_j, Te_1 = te_i, Tl_1 = tl_i,
            #                                          Te_2 = te_j, Tl_2 = tl_j, stringsAsFactors = FALSE))  # Commented out this block because rbind.data.frame becomes very slow when n is large.
        }
    }

    #return(list(Tm = Tm, Days = Days))
    return(Tm)
}

`%-%` <- function(t1, t2) as.integer(round(difftime(time1 = t1, time2 = t2, units = "days"), digits = 0))

.maxDays <- function(te_i, tl_i, te_j, tl_j, p_i, p_j) {
    # te_i, te_j: earliest possible collection date of isolates I and J, respectively
    # tl_i, tl_j: latest possible collection date of isolates I and J, respectively
    # Return: dt: day range of J - day range of I
    if (p_i == p_j && p_i == "YMD") {  # Explicit dates: te_i = tl_i, te_j = tl_j
        tm <- abs(tl_j %-% tl_i)
    } else {  # Calculate the maximum days from date ranges
        dl <- tl_j %-% tl_i  # dl = tl_j - tl_i
        de <- te_j %-% te_i  # de = te_j - te_i
        if (dl >= 0) {  # tl_j >= tl_i
            if (de >= 0) {  # tl_j > te_j >= te_i
                tm <- tl_j %-% te_i
            } else {  # te_j < te_i < tl_i <= tl_j
                tm <- max(c(tl_j %-% te_i, tl_i %-% te_j))
            }
        } else {  # tl_j < tl_i
            if (de >= 0) {  # tl_i > tl_j > te_j >= te_i
                tm <- max(c(tl_j %-% te_i, tl_i %-% te_j))
            } else {  # te_i > te_j
                tm <- tl_i %-% te_j
            }
        }
    }

    return(tm)
}
