#' @title Summarising phylogenetic distances for each clade in a rooted bifurcating phylogenetic tree
#' @description This function extracts clades from a rooted bifurcating phylogenetic tree and reports
#' the minimum, median, and the maximum of phylogenetic distances of each clade. Phylogenetic distances
#' may be SNP counts (neighbour-joining tree), number of substitutions per site (maximum-likelihood tree),
#' and so forth. The input tree will be midpoint rooted if it is unrooted.
#' @param tr An input bifurcating phylogenetic tree
#' @param dec Number of decimals kept for distance summaries
#' @returns A data frame summarising phylogenetic distances for each clade
#' @export
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 19 Feb 2022

distsPerClade <- function(tr, dec = 8) {
    require(ape)
    require(phytools)
    if (! is.rooted(tr)) {
        print("The input tree is unrooted and hence is rerooted using the midpoint.")
        tr <- ladderize(phy = midpoint.root(tr), right = FALSE)
    }
    taxa_num <- length(tr$tip.label)
    nodes_in <- (taxa_num + 1) : (taxa_num + tr$Nnode)  # IDs of internal nodes
    D <- do.call(rbind.data.frame, lapply(nodes_in, .summariseDistsPerClade, tr = tr))
    D$D_min <- round(D$D_min, digits = dec)
    D$D_max <- round(D$D_max, digits = dec)
    D$D_median <- round(D$D_median, digits = dec)
    return(D)
}

.summariseDistsPerClade <- function(i, tr) {  # i: index of an internal node; t: a phylo object
    c <- ape::extract.clade(phy = tr, node = i)  # The clade is a sub-tree rooted on the chosen internal node.
    return(data.frame(Node_in = i, Taxa_num = length(c$tip.label), D_min = min(c$edge.length),
                      D_max = max(c$edge.length), D_median = median(c$edge.length)))
}
