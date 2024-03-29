#' @title Plotting 95-percent confidence interval and dated tree in the result of function bactdate from the R package Bactdating.
#' @description Creating a tree figure with 95-percent confidence intervals of internal nodes' estimated dates. This function was modified
#' from function plot.resBactDating in github.com/xavierdidelot/BactDating/blob/master/R/methods.R.
#' @param x A dated tree generated by function bactdate of package Bactdating (github.com/xavierdidelot/BactDating).
#' @param years A sequence of years (integers) for ticks on the X axis. This vector can be created using the seq function.
#' @return An object can be plotted in R.
#' @export
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 3 Sep 2022; last update: 5 Sep 2022

plotCI <- function(x, years, ...) {  # Modified from plot.resBactDating in methods.R of the package
    require(ape)

    CI_range <- c(min(x$CI), max(x$CI))
    #xl <- range(pretty(CI_range - x$tree$root.time))
    xl <- round(range(years) - x$tree$root.time)  # Coordinates in the plot are based on the root time whose x is set to zero.

    # Draw the tree
    plot.phylo(x$tree, x.lim = xl, ...)  # This function does not draw any axis.
    #pre <- pretty(CI_range)
    #axis(side = 1, at = pre - x$tree$root.time, labels = pre)  # Re-label the X axis
    axis(side = 1, at = years - x$tree$root.time, labels = years)

    # Add tracks for CIs
    obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    transblue <- grDevices::rgb(0, 0, 1, 0.4)
    transred <- grDevices::rgb(1, 0, 0, 0.4)
    for(i in 1 : (Nnode(x$tree) + Ntip(x$tree)))
        if (x$CI[i, 1] != x$CI[i, 2]) {
            lines(x = c(x$CI[i, 1], x$CI[i, 2]) - x$tree$root.time,
                  y = rep(obj$yy[i], 2), lwd = 4,lend = 0,
                  col = ifelse(i <= Ntip(x$tree), transred, transblue))
        }
}
