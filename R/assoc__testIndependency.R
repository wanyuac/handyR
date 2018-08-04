#' @title Testing the independency of two variables according to a contingency table
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 21 Oct 2015

testIndependency <- function(m, method = "chisq", alt = "two.sided") {
    # m is an integer matrix representing a contingency table filled by columns
    # method = c("chisq", "fisher")
    # alt: is not applicable for chi-square tests

    if (method == "chisq") {
        ht <- chisq.test(m)
        d <- data.frame(chisq = ht$statistic, df = ht$parameter, p = ht$p.value)
    } else {    # alternatively, applies two-sided Fisher's exact tests
        ht <- fisher.test(x = m, conf.int = FALSE, alternative = alt)
        d <- data.frame(OR = ht$estimate, p = ht$p.value)
    }
    return(d)
}