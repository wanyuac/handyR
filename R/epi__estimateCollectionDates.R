#' @title Estimating dates of sample collection
#' @description This function estimates the range of collection dates "YYYY-MM" for each bacterial isolate based on a four-column data frame.
#' @param dates A data frame of four columns: Isolate, Collection_date (YYYY-MM-DD, YYYY, YYYY-MM), Collection_date_precision (Y, YM, YMD), Receipt_date (YYYY-MM-DD).
#' Unknown dates must be NA.
#' @return A data frame of the minimum and maximum collection dates of each isolate.
#' @export
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 17 August 2022; last update: 21 August 2022

estimateCollectionDates <- function(dates) {
    dates <- subset(dates, Collection_date_precision %in% c("YMD", "YM", "Y"))
    if (nrow(dates) > 0) {
        date_ranges <- do.call(rbind.data.frame, mapply(.dateRangesPerIsolate, dates$Isolate, dates$Collection_date, dates$Collection_date_precision, dates$Receipt_date,
                                                        MoreArgs = list(MONTH_DAYS = list("common year" = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
                                                                                          "leap year" = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))),
                                                        SIMPLIFY = FALSE, USE.NAMES = FALSE))
    } else {
        print("Error: Collection_date_format does not have any values of YMD, YM, or Y.")
        date_ranges <- NA
    }

    return(date_ranges)
}

.dateRangesPerIsolate <- function(i, c, c_type, r, MONTH_DAYS) {  # i: Isolate; c: collection date; c_type: collection-date type; r: receipt date, the date when the sample arrived at the lab.
    require(lubridate)
    use_receipt_date <- ! is.na(r)  # Otherwise, "if (days < 0)" and "if (as.integer(difftime(time1 = d_max, time2 = r)) > 0)" return an error of missing a TRUE/FALSE value.
    if (c_type == "YMD") {  # Nothing to do with a precise date
        d_min <- d_max <- c
    } else if (c_type == "YM") {
        yr <- as.integer(strsplit(x = c, split = "-", fixed = TRUE)[[1]][1])  # Year from "YYYY-MM"
        mon <- month(parse_date_time(x = c, orders = "ym"))  # mon (Month) is an integer.
        if (leap_year(yr)) {
            month_len <- MONTH_DAYS[["leap year"]][mon]
        } else {
            month_len <- MONTH_DAYS[["common year"]][mon]
        }
        d_min <- paste(c, "01", sep = "-")  # YYYY-MM-01
        d_max <- paste(c, as.character(month_len), sep = "-")
        if (use_receipt_date) {
            if (as.integer(difftime(time1 = d_max, time2 = r, units = "days")) > 0) {  # The isolate was received by the lab earlier than the estimated latest day. Class of function difftime's output: "difftime".
                d_max <- r
            }
        }
    } else {  # c_type = "Y"
        d_min <- paste(c, "01", "01", sep = "-")
        d_max <- paste(c, "12", "31", sep = "-")
        if (use_receipt_date) {
            receipt_year <- year(as.Date(r, format = "%Y-%m-%d"))
            collection_year <- as.integer(c)
            if (collection_year == receipt_year) {
                d_max <- r
            } else if (collection_year > receipt_year) {  # An impossible scenario indicating errors in dates
                print(paste("Error: receipt year of isolate", i, "was earlier than the collection year. Reset the date range to NA.", sep = " "))
                d_min <- d_max <- NA
            }
        }
    }

    # Validity check ###############
    if (use_receipt_date && (! is.na(d_min))) {
        days <- as.integer(difftime(time1 = r, time2 = d_min, units = "days"))  # days = time1 - time2. Normally, d_min always <= r.
        if (days < 0) {  # Suggests an error in dates: the isolate was received before its possible earliest collection date.
            print(paste("Error: receipt date of isolate", i, "was earlier than the earliest possible collection date. Reset the date range to NA.", sep = " "))
            d_min <- d_max <- NA
        }
    }

    return(data.frame(Isolate = i, Collection_date = c, Collection_date_precision = c_type, Receipt_date = r,
                      Collection_date_min = d_min, Collection_date_max = d_max, stringsAsFactors = FALSE))
}
