#' @title Table to data frame
#' @description Converts a table from function "table" to a two column data frame.
#' @param t Input table
#' @param new_names A character vector of two elements for new column names. Default: Name, Count.
#' @export
# (C) Copyright 2020 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Publication: 24 July 2020

tab2df <- function(t, new_names = NA) {
    df <- data.frame(Name = names(t), Count = as.integer(t), stringsAsFactors = FALSE)
    if (!is.na(new_names)[1]) {
        if (length(new_names) >= 2) {
            names(df) <- as.character(new_names)
        }
    }

    return(df)
}
