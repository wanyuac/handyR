#' @title Estimating dates of sample collection
#' @description This function estimates the range of collection dates "YYYY-MM" for each bacterial isolate based on a five-column data frame.
#' @param dates A data frame of five columns: Isolate, Collection_date (YYYY-MM-DD, YYYY, YYYY-MM), Collection_date_precision (Y, YM, YMD),
#' Receipt_date (YYYY-MM-DD), Date_of_birth (YYYY-MM-DD; set to NA if this information is not applicable or available). Unknown dates must be NA.
#' @return A data frame of the minimum and maximum collection dates of each isolate.
#' @export
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 17 August 2022; last update: 7 Sep 2022

estimateCollectionDates <- function(dates) {
    names(dates) <- c("Isolate", "Collection_date", "Collection_date_precision", "Receipt_date", "Date_of_birth")
    dates <- subset(dates, Collection_date_precision %in% c("YMD", "YM", "Y"))
    if (nrow(dates) > 0) {
        date_ranges <- do.call(rbind.data.frame, mapply(.dateRangesPerIsolate, dates$Isolate, dates$Collection_date, dates$Collection_date_precision, dates$Receipt_date, dates$Date_of_birth,
                                                        MoreArgs = list(MONTH_DAYS = list("common year" = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
                                                                                          "leap year" = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))),
                                                        SIMPLIFY = FALSE, USE.NAMES = FALSE))
    } else {
        print("Error: Collection_date_format does not have any values of YMD, YM, or Y.")
        date_ranges <- NA
    }

    return(date_ranges)
}

.dateRangesPerIsolate <- function(i, c, c_type, r, dob, MONTH_DAYS) {  # i: Isolate; c: collection date; c_type: collection-date type; r: receipt date, the date when the sample arrived at the lab; dob: date of birth
    # Estimates the collection date of isolate i.
    require(lubridate)
    use_receipt_date <- !is.na(r)  # Otherwise, "if (days < 0)" and "if (as.integer(difftime(time1 = d_max, time2 = r)) > 0)" return an error of missing a TRUE/FALSE value.
    use_DOB <- !is.na(dob)  # For each clinical isolate, the DOB sets the lower bound of the collection date.
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
        if (use_receipt_date) {  # The receipt date is the upper bound of the collection date.
            if (as.integer(difftime(time1 = d_max, time2 = r, units = "days")) > 0) {  # The isolate was received by the lab earlier than the estimated latest day. Class of function difftime's output: "difftime".
                d_max <- r
            }
        }
        if (use_DOB) {
            if (as.integer(difftime(time1 = dob, time2 = d_min, units = "days")) > 0) {  # If the date of birth is later than the current lower bound of possible collection dates
                d_min <- dob
            }
        }
    } else {  # c_type = "Y"
        d_min <- paste(c, "01", "01", sep = "-")
        d_max <- paste(c, "12", "31", sep = "-")
        collection_year <- as.integer(c)
        if (use_receipt_date) {
            receipt_year <- year(as.Date(r, format = "%Y-%m-%d"))
            if (collection_year == receipt_year) {
                d_max <- r
            } else if (collection_year > receipt_year) {  # An impossible scenario indicating errors in dates
                #print(paste("Error: receipt year of isolate", i, "was earlier than the collection year. Reset the date range to NA.", sep = " "))
                d_min <- d_max <- NA
            }
        }
        if (use_DOB) {
            birth_year <- year(as.Date(dob, format = "%Y-%m-%d"))
            if (collection_year == birth_year) {
                d_min <- dob
            } else if (collection_year < birth_year) {
                d_min <- d_max <- NA  # The metadata of this isolate may be wrong because the isolate cannot be collected before the birth of the patient.
            }
        }
    }

    # Verify metadata
    if (use_receipt_date && !is.na(d_min)) {
        if (as.integer(difftime(time1 = r, time2 = d_min, units = "days")) < 0) {  # days = time1 - time2. Normally, d_min always <= r.
            #print(paste("Error: receipt date of isolate", i, "was earlier than the earliest possible collection date. Reset the date range to NA.", sep = " "))
            d_min <- d_max <- NA
        }
    }
    if (use_DOB && !is.na(d_max)) {
        if (as.integer(difftime(time1 = d_max, time2 = dob, units = "days")) < 0) {  # Normally, d_max >= DOB. Isolates cannot be collected from a patient before their birth.
            d_min <- d_max <- NA  # Suggest users to check their metadata.
        }
    }

    return(data.frame(Isolate = i, Collection_date = c, Collection_date_precision = c_type, Date_of_birth = dob, Receipt_date = r,
                      Collection_date_min = d_min, Collection_date_max = d_max, stringsAsFactors = FALSE))
}
