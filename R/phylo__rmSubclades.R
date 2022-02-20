#' @title Removing IDs of sub-clades from a vector of clade names (internal nodes) of a bifurcating tree.
#' @description This function removes subclades from a list of clades and keeps non-overlapping clades.
#' This function is useful for processing the output of function distsPerClade of this package.
#' @param tr An input rooted bifurcating phylogenetic tree
#' @param nodes A character vector of internal node numbers or labels. Note that the characters should
#' be converted to integers if node numbers are used as input: as.integer(rmSubclades(...)).
#' @returns A vector of clades that do not overlap
#' @export
# (C) Copyright 2022 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 20 Feb 2022

rmSubclades <- function(tr, nodes) {
    require(ape)
    clade_tips <- list()
    clade_sizes <- integer()
    nodes_char <- character()
    nodes_pass <- character()
    use_node_numbers <- is.integer(nodes)

    # Summarise subtrees ===============
    for (i in nodes) {
        subtr <- extract.clade(phy = tr, node = i)
        if (use_node_numbers) {
            i <- as.character(i)
        }
        nodes_char <- append(nodes_char, i)
        clade_tips[[i]] <- subtr$tip.label
        clade_sizes[[i]] <- length(subtr$tip.label)
    }

    # Determining subclades ===============
    nodes_remain <- nodes_char
    for (i in nodes_char) {
        tips_i <- clade_tips[[i]]
        size_i <- clade_sizes[[i]]
        nodes_others <- nodes_remain[which(nodes_remain != i)]
        is_subclade <- FALSE  # Default status of tree i
        for (j in nodes_others) {
            tips_j <- clade_tips[[j]]
            size_j <- clade_sizes[[j]]
            if (size_i < size_j && setequal(intersect(x= tips_i, y = tips_j), tips_i)) {  # Set tips_i is a subset of tips_j.
                is_subclade <- TRUE  # Clade i is a subclade of j.
                break
            }
        }
        if (is_subclade) {
            nodes_remain <- nodes_others  # There's no need to compare a new clade to Clade i even if Clade i may contain subclades.
        } else {
            nodes_pass <- append(nodes_pass, i)
        }
    }

    # Return the result ===============
    if (use_node_numbers) {
        nodes_pass <- as.integer(nodes_pass)
    }
    return(nodes_pass)
}
