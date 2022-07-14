#' @title Summarising phylogenetic distances for each clade in a rooted bifurcating phylogenetic tree
#' @description This function extracts clades from a rooted bifurcating phylogenetic tree and reports
#' the minimum, median, and the maximum of phylogenetic distances of each clade. Phylogenetic distances
#' may be SNP counts (neighbour-joining tree), number of substitutions per site (maximum-likelihood tree),
#' and so forth. The input tree will be midpoint rooted if it is unrooted.
#' @param tr An input bifurcating phylogenetic tree
#' @param m (Optional) A distance matrix with row and column names matching to tip labels of the input tree
#' @param dec (Optional) Number of decimals kept for distance summaries (default: 8 digits)
#' @param boots (Optional) A data frame of bootstrap values. Format: Node (integer/character values), Bootstrap (floating-point numbers)
#' @returns A data frame summarising phylogenetic distances for each clade
#' @export
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 19 Feb 2022; last update: 13 July 2022

distsPerClade <- function(tr, m = NULL, dec = 8, boots = NULL) {
    require(ape)
    require(phytools)

    if (! is.rooted(tr)) {
        print("The input tree is unrooted and hence is rerooted using the midpoint.")
        tr <- ladderize(phy = midpoint.root(tr), right = FALSE)
    } else {
        print("The input tree is rooted. No tree manipulation will be performed.")
    }

    if ("node.label" %in% names(tr)) {
        nodes_in <- tr$node.label
    } else {
        taxa_num <- length(tr$tip.label)
        nodes_in <- (taxa_num + 1) : (taxa_num + tr$Nnode)  # IDs of internal nodes
    }

    D <- do.call(rbind.data.frame, lapply(nodes_in, .summariseDistsPerClade, tr = tr, m = m))
    D$D_min <- round(D$D_min, digits = dec)
    D$D_max <- round(D$D_max, digits = dec)

    if (! is.null(boots)) {
        names(boots) <- c("Node_in", "Bootstrap")
        D$Bootstrap <- boots$Bootstrap[match(D$Node_in, boots$Node_in)]
        D <- D[, c("Node_in", "Taxa_num", "Bootstrap", "D_min", "D_max")]
    }

    return(D)
}

.summariseDistsPerClade <- function(i, tr, m) {  # i: index of an internal node; tr: a phylo object; m: a distance matrix
    c <- ape::extract.clade(phy = tr, node = i)  # The clade is a sub-tree rooted on the chosen internal node.
    samples <- c$tip.label

    if (is.null(m)) {
        m_c <- ape::cophenetic.phylo(c)
    } else {
        m_c <- m[samples, samples]  # Extract distances from user's distance matrix
    }
    diag(m_c) <- NA  # Mask diagnoal entries

    return(data.frame(Node_in = i, Taxa_num = length(samples), D_min = min(m_c, na.rm = TRUE),
                      D_max = max(m_c, na.rm = TRUE)))
}
