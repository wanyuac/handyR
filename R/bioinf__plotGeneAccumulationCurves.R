#' @title Drawing gene-accumulation curves for a pangenome
#' @description This function is used for visualising results from function countPanGenes
#' @param L Result list of function countPanGenes
#' @param point_size Diameter of dots in the plot
#' @param line_width Width of curves
#' @param title Title of the plot (Default: "Gene accumulation")
#' @param axis_title_size Size of axis titles
#' @param axis_text_size Text size for axis ticks
#' @returns A ggplot2 object for plotting
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
# (C) Copyright 2023 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 15 Feb 2023; last update: 18 Feb 2023.

plotGeneAccumulationCurves <- function(L, point_size = 0.5, line_width = 0.5, title = "Gene accumulation",
                                       axis_title_size = 10, axis_text_size = 10) {
    require(ggplot2)
    n <- nrow(L$num$all)  # Number of isolates
    D <- rbind.data.frame(.MtoDf(L$num$hard_core, "Hard core"), .MtoDf(L$num$core, "Core"), .MtoDf(L$num$all, "Total"), stringsAsFactors = TRUE)
    avg_hardcore <- data.frame(x = 1 : n, y = L$avg$hard_core)
    avg_core <- data.frame(x = 1 : n, y = L$avg$core)
    avg_all <- data.frame(x = 1 : n, y = L$avg$all)
    p <- ggplot(data = D, mapping = aes(x = Genome, y = Gene, colour = Group)) +
        geom_point(size = point_size) +
        geom_line(data = avg_hardcore, mapping = aes(x = x, y = y), colour = "red", linewidth = line_width) +
        geom_line(data = avg_core, mapping = aes(x = x, y = y), colour = "orange", linewidth = line_width) +
        geom_line(data = avg_all, mapping = aes(x = x, y = y), colour = "blue", linewidth = line_width) +
        labs(x = "Number of genomes", y = "Number of genes", title = title) +
        scale_colour_manual(values = c("Hard core" = "#F88BC2", "Core" = "#FFE200", "Total" = "#47A9FA")) +
        theme_bw() +
        theme(axis.text = element_text(colour = "black", size = axis_text_size),
              axis.title = element_text(colour = "black", size = axis_title_size, face = "bold"),
              plot.title = element_text(colour = "black", size = axis_title_size, face = "bold"),
              legend.position = "none")
    return(p)
}

.MtoDf <- function(M, group = "Core") {
    # Converts a gene-count matrix to a data frame
    d <- data.frame(Genome = integer(0), Gene = integer(0), Group = character(0), stringsAsFactors = FALSE)
    perm <- ncol(M)
    for (x in 1 : nrow(M)) {
        d <- rbind.data.frame(d, data.frame(Genome = rep(x, times = perm),
                                            Gene = as.integer(M[x, ]),
                                            Group = rep(group, times = perm), stringsAsFactors = FALSE),
                              stringsAsFactors = FALSE)
    }
    return(d)
}
