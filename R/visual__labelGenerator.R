#' @title Creating labels for an axis at a given interval so that ticks can be denser than labels
#' @description For example, you may want the axis to be labelled by 0, 1k, 2k, etc.
#' @param ticks The ticks on an axis
#' @param interval N, a label is laid every N ticks
#' @param labels A vector of labels to be placed on the axis
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Last update: 29/4/2016

labelGenerator <- function(ticks, interval, labels){
    n.ticks <- length(ticks)
    newlabel <- character(n.ticks)  # a vector of empty characters
    j <- 1
    for (i in 1 : n.ticks) {
        if ((ticks[i] %% interval) == 0) {
            # a label should be laid at this tick
            newlabel[i] <- labels[j]
            j <- j + 1
        }
    }
    return(newlabel)
}