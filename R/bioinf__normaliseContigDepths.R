#' @title Median-based normalisation of contig depths
#' @description This function normalises fold coverages of contigs in a genome assembly. Depths that are close to the
#' median depth will be close to 1 after the normalisation.
#' @param tsv A Tab-delimited output (a TSV file) of four columns: isolate name, contig name, contig length,
#' and contig depth. This file can be generated from Shovill assemblies using script extract_contig_stats.sh in
#' package Assembly_toolkit (https://github.com/wanyuac/Assembly_toolkit).
#' @return A data frame with median-normalised depths
#' @author Yu Wan <wanyuac@@gmail.com>
#' @export
# (C) Copyright 2025 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# Release: 30 July 2025; last update: 30 July 2025.

normaliseContigDepths <- function(tsv) {
    require(tidyverse)
    contig_stats <- read_tsv(tsv, col_names = c("Isolate", "Contig", "Length", "Depth"), progress = FALSE, show_col_types = FALSE)
    median_depth <- median(contig_stats$Depth)
    IQR_depth <- IQR(contig_stats$Depth)
    result <- contig_stats %>% mutate(Depth_norm = round((Depth - median_depth) / IQR_depth + 1, digits = 1))
    return(result)
}
