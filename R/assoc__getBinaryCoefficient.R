#' @title Calculate association or similarity coefficient for two binary variables
#' and test for their independence.
#' @description A Jaccard similarity or phi association coefficient is computed
#' for binary variables a and b, such as a single pair of operational units (OUs).
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 21 Oct 2015

getBinaryCoefficient <- function(x, y, type = "Jaccard", hypo.test = "chisq", alt = "two.sided") {
    # calculates the coefficient for a pair of OUs as well as the probability of observing this contingency table
    # X and y must be two binary vectors of the same length
    # type = c("Jaccard", "phi", "p+J"), "p+J" means to return both the phi coefficients and Jaccard indices.
    # hypo.test: method = c("chisq", "fisher")

    # calculate the contingency table
    n <- length(x)
    if (n != length(y)) {
        r <- NA  # wrong inputs
    } else {
        a <- sum(x & y)  # e.g. c(0,1,0) & c(1,1,0) = c(FALSE, TRUE, FALSE)
        b <- sum(x & !y)
        c <- sum(!x & y)
        d <- n - a - b - c
        t <- matrix(data = c(a, c, b, d), nrow = 2, ncol = 2)  # defines the 2-by-2 contingency table
        p <- testIndependency(t, hypo.test, alt)$p
        if (type == "Jaccard") {
            r <- data.frame(Jc = a / (a + b + c), p = p, a = a, b = b, c = c, d = d)
        } else if (type == "phi") {
            r <- data.frame(phi = phiCoe(t), p = p, a = a, b = b, c = c, d = d)
        } else if (type == "p+J") {
            r <- data.frame(phi = phiCoe(t), Jc = a / (a + b + c), p = p, a = a, b = b, c = c, d = d)
        } else {  # More methods may be added in the future.
            r <- NA
        }
    }
    return(r)  # returns the coefficient and elements in the contingency table
}