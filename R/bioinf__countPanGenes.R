#' @title Counting genes in a pangenome when randomly sampling a series of genomes
#' @description This function is useful for creating gene-accumulation curves, also known as rarefaction curves.
#' @param P A binary presence-absence matrix with the dimension of n(genes) x n(isolates)
#' @param core Threshold defining hard/soft-core genes (Default: >=0.95)
#' @param hard_core Threshold defining hard-core genes (Default: >=0.99)
#' @param permutations Number of sampling for each sample size (Default: 20; minumum: 3)
#' @returns A list of two lists ("num" and "avg"): one for gene counts and the other for mean gene counts.
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
# (C) Copyright 2023 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# Release: 15 Feb 2023; last update: 18 Feb 2023.

countPanGenes <- function(P, core = 0.95, hard_core = 0.99, permutations = 20) {
    if (core > hard_core) {
        print("Warning: swapped parameters core and hard_core because the latter is smaller than the former")
        tmp <- core
        core <- hard_core
        hard_core <- tmp
    }
    if (permutations < 3) {
        print("Warning: reset parameter permutations to its minimum of three")
        permutations <- 3
    }
    num_isolates <- ncol(P)  # Number of isolates; assuming num_isolates >> 1
    indices <- 1 : num_isolates  # Indices of isolates
    iterations <- 1 : permutations  # Number of random sampling of isolates for estimating the variation in gene counts
    H <- matrix(0, nrow = num_isolates, ncol = permutations)  # Count matrix of hard-core genes; sample size per permutation x gene counts
    C <- matrix(0, nrow = num_isolates, ncol = permutations)  # Count matrix of hard/soft-core genes; sample size per permutation x gene counts
    A <- matrix(0, nrow = num_isolates, ncol = permutations)  # Count matrix of all genes; sample size per permutation x gene counts

    # Taking a single isolate (n = 1) at random, how many genes (0 or 1) are found?
    vH <- integer(0)  # A vector of hard-core gene counts
    vC <- integer(0)  # A vector of hard/soft-core gene counts
    vA <- integer(0)  # A vector of all-gene counts
    for (i in iterations) {
        j <- sample.int(n = num_isolates, size = 1, replace = FALSE)
        M <- as.integer(P[, j])  # In this iteration, M is a binary vector of a length of n(genes)
        g <- sum(M)  # Number of genes in this chosen isolate
        vH <- append(vH, g)  # Since only a single isolate was sampled, all genes are considered as hard-core genes in this sample.
        vC <- append(vC, g)  # Same as above
        vA <- append(vA, g)
    }
    H[1, ] <- vH
    C[1, ] <- vC
    A[1, ] <- vA

    # Taking 2, .., (num_isolates - 1) of isolates at random in each iteration
    for (n in 2 : (num_isolates - 1)) {
        vH <- integer(0)  # A vector of hard-core gene counts
        vC <- integer(0)  # A vector of hard/soft-core gene counts
        vA <- integer(0)  # A vector of all-gene counts
        for (i in iterations) {
            j <- sample.int(n = num_isolates, size = n, replace = FALSE)
            M <- P[, j]  # M is a matrix when j is a vector of a length > 1.
            gene_counts <- rowSums(M)  # Genetic occurrencies in these n genomes
            f <- gene_counts / n  # Genetic frequencies given the current sample of n isolates
            vH <- append(vH, sum(f >= hard_core))
            vC <- append(vC, sum(f >= core))
            vA <- append(vA, sum(gene_counts > 0))
        }
        H[n, ] <- vH
        C[n, ] <- vC
        A[n, ] <- vA
    }

    # The final iteration: taking all isolates (n = num_isolates)
    gene_counts <- rowSums(P)  # Surely, no variation will be seen in the gene content in sampled isolates
    f <- gene_counts / num_isolates
    H[num_isolates, ] <- rep(sum(f >= hard_core), times = permutations)
    C[num_isolates, ] <- rep(sum(f >= core), times = permutations)
    A[num_isolates, ] <- rep(sum(gene_counts > 0), times = permutations)  # Normally, sum(gene_counts) = nrow(P), but it is safer to use the sum function than using the latter.
    out <- list("num" = list("hard_core" = H, "core" = C, "all" = A),
                "avg" = list("hard_core" = round(rowMeans(H), digits = 2),
                             "core" = round(rowMeans(C), digits = 2),
                             "all" = round(rowMeans(A), digits = 2)))
    return(out)
}
